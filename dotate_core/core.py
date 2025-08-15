# Dotate - A tool for annotating ECOD protein domains based on HMM profile search
# Copyright (C) 2025 ELSI - Andre Lecona
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

import os
import json
import sqlite3
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool
from sqlalchemy import create_engine


def clean_dataframe(df, hmmN, i_eN):
    """Cleans and filters the input DataFrame based on i_Evalue and HMM coverage.

    Args:
        df (pd.DataFrame): Input DataFrame containing domain information.
        hmmN (float): HMM coverage threshold.
        i_eN (float): i_Evalue threshold.

    Returns:
        pd.DataFrame: Filtered and sorted DataFrame with relevant columns.
    """

    columns_to_convert = ['i_Evalue', 'from_env', 'to_env', 'qlen', 'from_hmm', 'to_hmm', 'tlen'] #convert to float()
    df.loc[:, columns_to_convert] = df[columns_to_convert].apply(pd.to_numeric, errors='coerce')
    df['hmm_cov'] = (df['to_hmm'] - df['from_hmm'] + 1) / df['qlen']
    filtered_df = df[(df['i_Evalue'] < i_eN) & (df['hmm_cov'] > hmmN)].copy()
    
    # Select only the desired columns and sort by 'from_env' in ascending order
    filtered_columns_df = filtered_df[['target_name', 'tlen', 'query_name', 'from_env', 'to_env', 'i_Evalue', 'hmm_cov']] \
        .sort_values(by='i_Evalue', ascending=True)
    
    return filtered_columns_df

def calculate_overlap_percentage(test_row, accepted_row):
    """
    Calculates the percentage of overlap between two domain ranges.

    Args:
        test_row (pd.Series): A row from the DataFrame representing the domain to test.
        accepted_row (pd.Series): A row from the accepted domains representing a previously processed domain.

    Returns:
        float: The overlap percentage of the test domain relative to its length.
    """
    # Extract the ranges
    start1, end1 = test_row['from_env'], test_row['to_env']
    start2, end2 = accepted_row['from_env'], accepted_row['to_env']
    
    # Find the overlap range
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    
    # Calculate overlap length and test_row domain length
    overlap_length = max(0, overlap_end - overlap_start + 1)
    length1 = end1 - start1 + 1
    
    # Return overlap percentage relative to test_row
    return overlap_length / length1 if length1 > 0 else 0

def add_gap_rows(df, amino):
    """
    Adds gap rows to the DataFrame to cover unannotated regions between existing domains.

    Args:
        df (pd.DataFrame): Input DataFrame containing the domain annotations.
        amino (int): Minimum gap size (in amino acids) required to insert a gap row.

    Returns:
        pd.DataFrame: The DataFrame with added gap rows to fill unannotated regions.
    """
    # Ensure the ranges are sorted by 'from_env'
    df = df.sort_values(by='from_env').reset_index(drop=True)
    new_rows = []

    def add_gap_row(start, end):
        return {
            'target_name': df.iloc[0]['target_name'],
            'tlen': 0,
            'query_name': 'UNN',
            'from_env': start,
            'to_env': end,
            'i_Evalue': 0,
            'hmm_cov': 0,
            'domain_cov': 0
        }

    # Check for gaps and create new rows
    if df.iloc[0]['from_env'] > amino:
        new_rows.append(add_gap_row(0, df.iloc[0]['from_env'] - 1))
    
    for i in range(1, len(df)):
        gap_size = df.iloc[i]['from_env'] - df.iloc[i-1]['to_env'] - 1
        if gap_size >= amino:
            new_rows.append(add_gap_row(df.iloc[i-1]['to_env'] + 1, df.iloc[i]['from_env'] - 1))
    
    if df.iloc[-1]['tlen'] - df.iloc[-1]['to_env'] >= amino:
        new_rows.append(add_gap_row(df.iloc[-1]['to_env'] + 1, df.iloc[-1]['tlen']))

    # Return the merged DataFrame
    return pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True).sort_values(by='from_env').reset_index(drop=True)

def filter_accepted_rows(df, coverage, amino):
    """
    Filters and processes the input DataFrame by removing overlapping domains and adding gap rows.

    Args:
        df (pd.DataFrame): Input DataFrame containing domain information.
        coverage (float): Minimum coverage threshold for domain overlap.
        amino (int): Minimum gap size (in amino acids) required to insert a gap row.

    Returns:
        pd.DataFrame: Filtered and sorted DataFrame with relevant columns, including gap rows where necessary.
    """
    accepted_rows = []

    for i, test_row in df.iterrows():
        # Store the original length of the domain
        initial_length = test_row['to_env'] - test_row['from_env'] + 1
        
        # Check for overlaps with accepted rows
        overlaps = [calculate_overlap_percentage(test_row, accepted_row) for accepted_row in accepted_rows]

        if all(overlap < 1 - coverage for overlap in overlaps):

            # Handle overlap reduction if necessary
            for accepted_row in accepted_rows:
                overlap_percentage = calculate_overlap_percentage(test_row, accepted_row)
                if overlap_percentage > 0:  # Overlap found, but less than threshold
                    # Adjust to eliminate overlap
                    start1, end1 = test_row['from_env'], test_row['to_env']
                    start2, end2 = accepted_row['from_env'], accepted_row['to_env']
                    if start1 < start2:
                        test_row['to_env'] = start2 - 1
                    else:
                        test_row['from_env'] = end2 + 1

            # Add domain coverage to the test_row
            test_row['domain_cov'] = (test_row['to_env'] - test_row['from_env'] + 1) / initial_length
            accepted_rows.append(test_row)

    # Add gap rows and finalize results
    result = add_gap_rows(pd.DataFrame(accepted_rows), amino).drop('tlen', axis=1)
    result['hmm_cov'] = result['hmm_cov'].round(3)
    result['domain_cov'] = result['domain_cov'].round(3)

    return result

def process_gene(df_gene, hmm_cov_co, iE_score_co, domain_cov_co, unanotated_co):
    """Processes a single gene by extracting, cleaning, and annotating its protein.

    Args:
        gene (str): Gene identifier.
        df_genes (dict): Dictionary of DataFrames, grouped by 'target_name'.
        hmm_cov_co (float): HMM coverage cutoff.
        iE_score_co (float): i-Evalue cutoff.
        domain_cov_co (float): Domain coverage cutoff.
        unanotated_co (int): Minimum length for unannotated regions.

    Returns:
        pd.DataFrame: Annotated protein DataFrame.
    """

    df_cleaned = clean_dataframe(df_gene, hmm_cov_co, iE_score_co)

    if df_cleaned.empty:
        tlen = df_gene['tlen'].iloc[0]
        unanotated_row = {
            'target_name': df_gene.iloc[0]['target_name'],
            'domain': 'UNN',
            'from_env': 0,
            'to_env': tlen,
            'i_Evalue': float('nan'),
            'hmm_cov': float('nan'),
            'domain_cov': float('nan'),
        }
        return pd.DataFrame([unanotated_row])
    
    return filter_accepted_rows(df_cleaned, domain_cov_co, unanotated_co).rename(columns={'query_name': 'domain'})

def ECODmapping(df):
    """
    Maps the 'domain' column in the input DataFrame to the corresponding ECOD IDs
    based on develop289 version, and renames the column to 'f_id'.

    Args:
        df (pandas.DataFrame): The input DataFrame containing a 'query_name' column.
        mapping_file (str): The path to the JSON file containing the ECOD mapping,

    Returns:
        pandas.DataFrame: The updated DataFrame 
    """
    #with open('nm2id.json', 'r') as file:
    #    f2x_mapping = json.load(file)
    
    from dotate_core.utils import f2x_mapping

    # Map the 'domain' column to corresponding ECOD IDs from the mapping data
    df['f_id'] = df['domain'].map(f2x_mapping).fillna('UNN')
    columns_order = ['target_name', 'f_id', 'domain', 'from_env', 'to_env', 'i_Evalue', 'hmm_cov', 'domain_cov']

    # Return the updated DataFrame
    return df[columns_order]

def worker(args):
    """
    Processes a single gene group by running the `process_gene` function.

    Args:
        args (tuple): Contains the following elements:
            - target_name (str): Name of the target gene.
            - gene_group (pd.DataFrame): DataFrame containing the gene's domain data.
            - hmm_cov_co (float): HMM coverage threshold.
            - iE_score_co (float): E-value threshold for domain inclusion.
            - domain_cov_co (float): Minimum domain coverage required.
            - unanotated_co (int): Minimum gap size for unannotated regions.

    Returns:
        pd.DataFrame: Processed gene group with filtered and annotated domain data.
    """
    target_name, gene_group, hmm_cov_co, iE_score_co, domain_cov_co, unanotated_co = args
    return process_gene(gene_group, hmm_cov_co, iE_score_co, domain_cov_co, unanotated_co)

def dotate(hmm_file, mapping=False, hmm_cov_co=0.75, iE_score_co=0.01, domain_cov_co=0.75, unanotated_co=10, cores=1, chunksize=100):
    """
    Reads an HMM search result file, processes gene groups in parallel, and returns a consolidated DataFrame.

    Args:
        hmm_file (str): Path to the HMM search results file.
        hmm_cov_co (float, optional): HMM coverage threshold. Default is 0.75.
        iE_score_co (float, optional): Inclusion E-value threshold. Default is 0.01.
        domain_cov_co (float, optional): Minimum domain coverage required. Default is 0.75.
        unanotated_co (int, optional): Minimum gap size (in amino acids) to insert a gap row. Default is 10.
        cpu (int, optional): Number of CPU cores to use for multiprocessing. Default is 1.
        chunksize (int, optional): Number of gene groups to process per worker in a batch. Default is 70.

    Returns:
        pd.DataFrame: A concatenated DataFrame of processed gene groups with domain annotations.
    """

    # Define column names for the HMM search file
    column_names = [
        'target_name', 'accession', 'tlen', 'query_name', 'query_accession', 'qlen',
        'E_value', 'score_full', 'bias_full', 'n', 'N', 'c-Evalue', 'i_Evalue', 'score_domain', 'bias_domain',
        'from_hmm', 'to_hmm', 'from_ali', 'to_ali', 'from_env', 'to_env']

    df = pd.read_csv(hmm_file, delimiter=r'\s+', comment='#', header=None, engine='python', usecols=range(21), names=column_names)

    # Group by 'target_name' and prepare input for multiprocessing
    df_genes = [(target_name, gene_group, hmm_cov_co, iE_score_co, domain_cov_co, unanotated_co) 
                for target_name, gene_group in df.groupby('target_name')]

    # Process gene groups in parallel using multiprocessing
    with Pool(processes=cores) as pool:
        proteome = list(tqdm(pool.imap_unordered(worker, df_genes, chunksize), total=len(df_genes)))

    # Concatenate the processed gene groups into a single DataFrame
    proteome_df = pd.concat(proteome, ignore_index=True)
    proteome_df.reset_index(drop=True, inplace=True)

    if mapping:
        return ECODmapping(proteome_df)
    else: return proteome_df

def store_SQL(hmm_file: str, df: pd.DataFrame, db_path: str):
    """
    Stores the given DataFrame into an SQLite database file. If the file does not exist, it is created.
    
    Parameters:
    - hmm_file (str): Name of the HMM file (used as the table name).
    - df (pd.DataFrame): DataFrame to store.
    - db_path (str): Path to the SQLite database file (e.g., 'my_db.sqlite').
    """
    if not db_path:
        raise ValueError("A valid SQLite database file path must be provided.")

    # Create SQLAlchemy engine for SQLite
    engine = create_engine(f"sqlite:///{db_path}")

    # Create a valid table name from the hmm file
    table_name = os.path.splitext(os.path.basename(hmm_file))[0]
    table_name = table_name.replace(".", "_").replace(" ", "_").lower()

    # Store the DataFrame in SQLite
    df.to_sql(name=table_name, con=engine, if_exists="replace", index=False)

    print(f"Data stored successfully in table '{table_name}' in SQLite database at '{db_path}'")


def store_fasta(df: pd.DataFrame, output_file: str):
    """
    Writes a FASTA file grouping domain IDs by protein.
    
    Parameters:
        df (pd.DataFrame): Input DataFrame with columns ['target_name', 'domain', 'from_env', 'to_env'].
        output_file (str): Path to the output FASTA file.
    """
    df = df.sort_values(by=['target_name', 'from_env'])

    if 'f_id' in df.columns:
        col='f_id'
    else:
        col='domain'

    def format_col(group):
        formatted_cols = []
        for _, row in group.iterrows():
            if row[col] == 'UNN':
                length = row['to_env'] - row['from_env'] + 1
                formatted_cols.append(f'({length})')
            else:
                formatted_cols.append(row[col])
        return '-'.join(formatted_cols)
    
    grouped = df.groupby('target_name', group_keys=False)[[col, 'from_env', 'to_env']].apply(format_col).reset_index(name=col)

    with open(output_file, 'w') as fasta:
        for i, row in grouped.iterrows():
            fasta.write(f">{row['target_name']}\n")
            fasta.write(f"{row[col]}\n")
            fasta.write("\n")

if __name__ == "__main__":
    FASTA = 'testdotate.fasta'
    HMM = 'testdotate.tbl'
    SQL = 'SQL.db'

    try:
        proteome_df = dotate(hmm_file=HMM)        
        unannotated_count = proteome_df['domain'].eq('UNN').sum()
    except Exception as e:
        print(f"Error during mapping: {e}")
        exit(1)


    try:
        store_fasta(ECODmapping(proteome_df), FASTA)
        print("FASTA complete")
    except Exception as e:
        print(f"Error during fasta creation: {e}")


    try:
        output_file = 'dotate_out.tsv'
        ECODmapping(proteome_df).to_csv(output_file, sep='\t', index=False, header=True)
        print("TSV complete")
    except Exception as e:
        print(f"Failed to save TSV file: {e}")


    try:
        store_SQL(HMM, ECODmapping(proteome_df), SQL)
        print("SQL complete")
    except Exception as e:
        print(f"Failed to store SQL: {e}")

    print("Process completed.")