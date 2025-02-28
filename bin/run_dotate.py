#!/usr/bin/env python

# Dotate - A tool for annotating ECOD protein domains based on HMM profile search
# Copyright (C) 2024 ELSI - Andre Lecona
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


import sys
import os
import argparse
import json
from dotate_core.core import dotate, store_SQL  # Import the main function and store_SQL

def load_sql_config(sql_config_path):
    """Load the SQL configuration from a JSON file."""
    with open(sql_config_path, 'r') as config_file:
        return json.load(config_file)

def main():
    print("Starting DOTATE process...")

    # Set up argument parsing
    parser = argparse.ArgumentParser(
        description="Dotate: A tool for annotating protein domains based on HMM domain table output."
    )
    
    # Command-line arguments
    parser.add_argument('hmm_file', type=str, help='Path to the HMM domain table output file')
    parser.add_argument('--ECODmapping', action='store_true', help='Enable mapping (default is False)')
    parser.add_argument('--hmm_cov_co', type=float, default=0.75, help='HMM coverage cutoff (default is 0.75)')
    parser.add_argument('--iEvalue_co', type=float, default=0.01, help='iE-value score cutoff (default is 0.01)')
    parser.add_argument('--domain_cov_co', type=float, default=0.75, help='Domain coverage cutoff (default is 0.75)')
    parser.add_argument('--unanotated_co', type=int, default=10, help='Unannotated coverage threshold (default is 10)')
    parser.add_argument('--cores', type=int, default=1, help='Number of cores to use (default is 1)')
    parser.add_argument('--chunksize', type=int, default=100, help='Chunk size for processing (default is 100)')
    parser.add_argument('--SQL', type=str, help='Path to SQLconfig.json file.')
    parser.add_argument('--version', action='version', version='Dotate 0.1.0 (GPLv3)')

    # Parse arguments
    args = parser.parse_args()

    # Print the help message if --help is called
    if '--help' in sys.argv:
        parser.print_help()
        return

    HMM = args.hmm_file
    print(f"Processing HMM file: {HMM}")

    try:
        # Call the main function with the parsed arguments
        proteome_df = dotate(
            hmm_file=HMM,
            mapping=args.ECODmapping,
            hmm_cov_co=args.hmm_cov_co,
            iE_score_co=args.iEvalue_co,
            domain_cov_co=args.domain_cov_co,
            unanotated_co=args.unanotated_co,
            cores=args.cores,
            chunksize=args.chunksize
        )

        # Count unannotated entries
        unannotated_count = proteome_df['f_id'].eq('Unannotated').sum()
        print(f"Annotation complete. Found {unannotated_count} 'Unannotated' entries.")

    except Exception as e:
        print(f"Error during mapping: {e}")
        exit(1)

    try:
        if args.SQL:
            # Load the SQLconfig from the provided JSON file
            try:
                SQLconfig = load_sql_config(args.SQL)
                print("SQL configuration found. Storing data...")
                store_SQL(HMM, proteome_df, SQLconfig)
                print(f"Data successfully stored in database '{SQLconfig.get('database', 'default_database')}'.")

            except Exception as e:
                print(f"Error loading SQL config: {e}")
                exit(1)

        else:
            # Save the output to a .dotate.tsv file
            output_file = f"{os.path.splitext(HMM)[0]}.dotate.tsv"
            proteome_df.to_csv(output_file, sep='\t', index=False, header=True)
            print(f"Results saved to '{output_file}'. Process completed.")

    except Exception as e:
        print(f"Unexpected error while storing data in SQL or saving TSV: {e}")
        exit(1)

    print("Process completed.")

if __name__ == "__main__":
    main()
