#!/usr/bin/env python

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
import argparse
import json
import multiprocessing
from dotate_core.core import dotate, store_SQL, store_fasta

def load_sql_config(sql_config_path):
    """Load the SQL configuration from a JSON file."""
    if not os.path.isfile(sql_config_path):
        raise argparse.ArgumentTypeError(f"SQL config file not found: {sql_config_path}")
    try:
        with open(sql_config_path, 'r') as config_file:
            return json.load(config_file)
    except json.JSONDecodeError:
        raise argparse.ArgumentTypeError(f"Invalid JSON format: {sql_config_path}")

def validate_positive(value):
    ivalue = float(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(f"Invalid value: {value}. Must be a positive number.")
    return ivalue

def validate_float(value):
    ivalue = float(value)
    if ivalue <= 0 and ivalue > 1:
        raise argparse.ArgumentTypeError(f"Invalid value: {value}. Must be a positive number between 0 and 1.")
    return ivalue

def validate_path(path):
    directory = os.path.dirname(path)
    if directory and not os.path.exists(directory):
        raise argparse.ArgumentTypeError(f"Invalid path: {path}. Directory does not exist.")
    return path

def validate_file(path):
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"File not found: {path}")
    return path

def get_num_cores(cores):
    if cores == -1:
        return max(1, multiprocessing.cpu_count() - 1)
    if cores < 1:
        raise argparse.ArgumentTypeError("Cores must be at least 1 or -1 for auto-selection.")
    return cores

def main():
    print("Starting DOTATE process...")

    # Set up argument parsing
    parser = argparse.ArgumentParser(
        description="Dotate: A tool for annotating protein domains based on HMM domain table output."
    )
    
    # Command-line arguments
    parser.add_argument('hmm_file', type=validate_file, help='Path to the HMM domain table output file')
    parser.add_argument('--ECODmapping', action='store_true', help='Enable mapping (default is False)')
    parser.add_argument('--hmm_cc', type=validate_float, default=0.75, help='HMM coverage cutoff (Values 0 to 1 -- default is 0.75)')
    parser.add_argument('--iEvalue_co', type=validate_positive, default=0.01, help='iE-value score cutoff (default is 0.01)')
    parser.add_argument('--domain_cc', type=validate_float, default=0.75, help='Domain coverage cutoff (Values 0 to 1 default is 0.75)')
    parser.add_argument('--unn_co', type=validate_positive, default=50, help='Length cutoff for detecting unannotated segments (default is 50)')
    parser.add_argument('--cores', type=int, default=1, help='Number of cores to use (-1 for max-1, default is 1)')
    parser.add_argument('--chunksize', type=int, default=100, help='Chunk size for processing (default is 100)')
    parser.add_argument('--sql', type=validate_file, help='Path to SQLconfig.json file.')
    parser.add_argument('--fasta', type=validate_path, help='Path to output FASTA file.')
    parser.add_argument('--tsv', type=validate_path, help='Path to output TSV file.')
    parser.add_argument('--version', action='version', version='Dotate 1.1.1 (GPLv3)')

    # Parse arguments
    args = parser.parse_args()
    
    # Adjust cores if -1 is provided
    args.cores = get_num_cores(args.cores)

    HMM = args.hmm_file
    print(f"Processing HMMsearch file: {HMM}")

    try:
        # Call the main function with the parsed arguments
        proteome_df = dotate(
            hmm_file=HMM,
            mapping=args.ECODmapping,
            hmm_cov_co=args.hmm_cc,
            iE_score_co=args.iEvalue_co,
            domain_cov_co=args.domain_cc,
            unanotated_co=args.unn_co,
            cores=args.cores,
            chunksize=args.chunksize
        )

        # Count unannotated entries
        unannotated_count = proteome_df['domain'].eq('UNN').sum()
        annotated_count = proteome_df['domain'].ne('UNN').sum()
        print(f"Annotation complete. Found {annotated_count} domains and {unannotated_count} 'Unannotated' entries")

    except Exception as e:
        print(f"Error during mapping: {e}")
        exit(1)

    try:
        if args.sql:
            try:
                SQLconfig = load_sql_config(args.sql)
                print("SQL configuration found. Storing data...")
                store_SQL(HMM, proteome_df, SQLconfig)
                print(f"Data successfully stored in database '{SQLconfig.get('database', 'default_database')}'.")
            except Exception as e:
                print(f"Error loading SQL config: {e}")

        # Save to specified output formats
        if args.tsv:
            try:
                proteome_df.to_csv(args.tsv, sep='\t', index=False, header=True)
                print(f"Results saved to '{args.tsv}'.")
            except Exception as e:
                print(f"Error saving TSV file: {e}")

        if args.fasta:
            try:
                store_fasta(proteome_df, args.fasta)
                print(f"FASTA file saved to '{args.fasta}'.")
            except Exception as e:
                print(f"Error saving FASTA file: {e}")
          
        elif not args.sql and not args.tsv and not args.fasta:
            try:
                default_tsv = f"{os.path.splitext(HMM)[0]}.dotate.tsv"
                proteome_df.to_csv(default_tsv, sep='\t', index=False, header=True)
                print(f"Results saved to '{default_tsv}'.")
            except Exception as e:
                print(f"Error saving default TSV file: {e}")

    except Exception as e:
        print(f"Unexpected error while storing data: {e}")
        exit(1)

    print("Process completed.")

if __name__ == "__main__":
    main()