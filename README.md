# Dotate: Annotate Protein Domains with HMM Output

Dotate is a tool for annotation of protein domains based on HMMsearch results. It is particularly useful as part of a bioinformatics proteomics pipeline, enabling the annotation of multiple proteomes or large protein databases. Dotate also supports integration with an SQL database for storing results. When working with ECOD HMM profiles, Dotate can use ECOD's F-ID notation for domain groups.

## Features:
- Annotate protein domains using from an HMMsearch domain-table out file.
- Stores results in an SQL database, a TSV file, and/or a FASTA-like file.
- Optimized for large datasets.

## Usage:
dotate <HMMSEARCH_OUTFILE> [--help] [--ECODmapping] [--hmm_cc HMM_CC] [--iEvalue_co IEVALUE_CO] [--domain_cc DOMAIN_CC] [--unn_co UNN_CO] [--cores CORES] [--chunksize CHUNKSIZE] [--sql SQL] [--fasta FASTA]

## Example:
(1) Independently run HMMER search.
- hmmsearch --domtblout result.tbl your_model.hmm your_sequences.fasta

(2) Use DOTATE to annotate and parse your result.  
Annotate a protein domain and store in SQL:
- dotate result.tbl --hmm_cc 0.8 --iEvalue_co 0.005 --SQL my_database.json --ECODmapping

Annotate and store as TSV file (default):
- dotate result.tbl --hmm_cc 0.8 --iEvalue_co 0.005 --ECODmapping

## SQLconfig
Example for your SQLconfig.json file:

```json
{
    "host": "000.000.00",
    "user": "user",
    "password": "password",
    "database": "database"
}
```

## Installation

Clone the repository and install the package:

```bash
git clone https://github.com/AndreLecona/Dotate.git
pip install ./dotate
```

## License
Dotate is licensed under the GNU General Public License v3.0 (GPLv3).  
See the [LICENSE](LICENSE) file for more details.