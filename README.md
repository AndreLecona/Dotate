# Dotate: Annotate Protein Domains with HMM Output

Dotate is a bioinformatics tool for annotating protein domains based on HMM domain table output. It is designed for searching ECOD domains and storing the results in a SQL database.

## Features:
- Annotate protein domains using HMM domain tables.
- Store results in a SQL database or as a TSV file.
- Optimized for large datasets with chunk processing.

## Usage:
dotate <HMMSEARCH_OUTFILE> [--help] [--ECODmapping] [--hmm_cov_co HMM_COV_CO] [--iEvalue_co IEVALUE_CO] [--domain_cov_co DOMAIN_COV_CO] [--unanotated_co UNANOTATED_CO] [--cores CORES] [--chunksize CHUNKSIZE] [--SQL SQL]

## Example:
Annotate a protein domain and store in SQL:
- dotate testdotate.tbl --hmm_cov_co 0.8 --iE_score_co 0.005 --SQL my_database.db --ECODmapping

Annotate and store as TSV file:
- dotate testdotate.tbl --hmm_cov_co 0.8 --iE_score_co 0.005 --ECODmapping

## SQLconfig
Example for your SQLconfi.json file:

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
cd dotate
pip install -e .
```

## License
Dotate is licensed under the GNU General Public License v3.0 (GPLv3).  
See the [LICENSE](LICENSE) file for more details.