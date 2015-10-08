#!/usr/bin/env python

'''
Author: Prag Batra

Purpose:
    
    Provides a single place to keep track of all constants used in STMP.

Explanation:

    See purpose above. Note: all column names should be lowercase
    for comparison purposes.

Example:

    For usage, see stmp2.py where these methods are all called.

'''

import collections
from collections import OrderedDict

# annotation_db_file = 'data/annotationDB.sqlite'
# db_beds_dir = 'data/db_beds'
tiering_db_file = 'data/tieringDB.sqlite'

# must exactly match the given column name
# NOTE: the standard VCF columns (chrom, start, id, ref, alt, qual, filter, info) are currently hardcoded so this dictinary is ignored for these values
col_table_types_absolute = { 
    'chr': 'varchar(10)',
    'chrom': 'varchar(10)',
    'start': 'int',
    'stop': 'int',
    'ref': 'varchar(512)',
    'alt': 'varchar(512)',
    'pop_freq_max': 'float',  
    'info': 'varchar(255)',
    'format': 'varchar(255)'
                   }

# Only needs to match the beginning of the given column name.
# These are listed in search order (first match is the one used).
col_table_types_prefixes = OrderedDict([
    ('shuffle_wellderly', 'varchar(255)')
                                        ])

# { 
#     'shuffle_wellderly': 'varchar(255)'
#                          
#                          }

