# Licensed under GPLv3 - see LICENSE

"""
Dotate: A tool for annotating ECOD protein domains based on HMM profile search

Provides functionality to search ECOD (Evolutionary Classification of Protein Domains) domains
using HMMsearch domain-table output files and stores the resulting annotations in an SQL database.

Main Features:
- Annotate protein domains based on HMM search results.
- Supports searching ECOD domains.
- Store domain annotations in a database for later querying.

Version: 1.1.1
Author: Andre Lecona Buttelli
Email: andrelecona@elsi.jp
"""

#metadata
__version__ = '1.1.1'
__author__ = 'Andre Lecona Buttelli'
__author_email__ = 'andrelecona@elsi.jp'
__description__ = 'Dotate: A bioinformatics tool for annotating protein domains based on HMMsearch domain-table output.'
__url__ = 'https://github.com/AndreLecona/Dotate'

# Expose important functionality for users
from .core import dotate
from .utils import f2x_mapping

__all__ = [
    'dotate',
    'f2x_mapping',
]