#!/usr/bin/env python

'''
Author: Prag Batra prag@stanford.edu

Purpose:
    Contains helper functions to run preflight checks prior to running the entire STMP annotation pipeline (e.g. integrity checks on column names).
    

Explanation:

    

Example:

    

'''

import db_utils
import yaml_keys
import yaml_utils
import stmp_annotation_util
import vcfHeaders
import re

# gets a list (dictionary) of all tiering columns in the yaml (for checking)
def get_all_tiering_cols(yaml_commands):
    modules_yaml = yaml_commands[yaml_keys.kModules]
    tiering_yaml = modules_yaml[yaml_keys.kTiering]
    tiering_cols = {}
    tiering_cols[tiering_yaml[yaml_keys.kTGeneNameCol]] = 1
    for col in tiering_yaml[yaml_keys.kTAlleleFreqCols]:
        if(col is list):
            for subcol in col:
                tiering_cols[subcol] = 1
    for col in tiering_yaml[yaml_keys.kTFunctionalCols]:
        if(col is list):
            for subcol in col:
                tiering_cols[subcol] = 1
    for col in tiering_yaml[yaml_keys.kTConservationCols]:
        if(col is list):
            for subcol in col:
                tiering_cols[subcol] = 1
    for col in tiering_yaml[yaml_keys.kTPathogenicityCols]:
        if(col is list):
            for subcol in col:
                tiering_cols[subcol] = 1
    for col in tiering_yaml[yaml_keys.kTToleranceZScoreCols]:
        if(col is list):
            for subcol in col:
                tiering_cols[subcol] = 1
    
    return tiering_cols

# checks column names against the database for both annotation and tiering. Identifies inconsistencies bw what is actually loaded into the db compared to the db spec (from datasets.yml).
def preflight_checks(yaml_commands, db_conn, raise_exception=True):
    errors = []
    warnings = []
    c = db_conn.cursor()
    modules_yaml = yaml_commands[yaml_keys.kModules]
    
    # TODO: 
    # - add more checks for any currently unchecked yaml keys, etc. in modules.yml and datasets.yml (e.g. throw an error if a yaml key like "Consensus_Columns" is missing) 
    # - if there are any currently unchecked dataset column references in modules.yml, check those (cross-reference against datasets.yml by using yaml_utils.convertColumns() to ensure the appropriate columns can be found in datasets.yml)
    # - add checks for all external files (e.g. target gene lists for tiering). Warn/error if target gene list files are empty/missing.
    # - add checks for snpeff and annovar path
    # - check that summary column (if used, e.g. as gene name column for tiering) is specified with appropriate suffix (same as the suffix specified as Consensus_Column_Suffix)
    # - perhaps even do a simulated annotation run (just the header) and compare absolute headers in modules.yml against that, too
    
    # WARNING: currently not checking absolute columns (only columns specified with dot (.) notation)
    # 1) check downstream column references in modules.yml against datasets.yml
    # 1-1) check consensus columns
    consensus_cols_yaml = modules_yaml[yaml_keys.kAnnotation][yaml_keys.kAConsensusColumns]
    for consensus_col in consensus_cols_yaml:
        if consensus_col != yaml_keys.kAConsensusColumnsOrder:
            yaml_utils.convertColumns(consensus_cols_yaml[consensus_col], yaml_commands) # will raise exception if can't find one of the columns with dot (.) notation
    # 1-2) check tiering columns
    tiering_cols_yaml = get_all_tiering_cols(yaml_commands)
    yaml_utils.convertColumns(tiering_cols_yaml, yaml_commands) # will raise exception if can't find one of the columns with dot (.) notation
    
    # 2) check datasets.yml against database (and vice-versa)
    # 2-1) check table names to see if any are missing from db/datasets.yml
    datasets_yaml = yaml_utils.get_datasets(yaml_commands)
#     dataset_annotation_names_dict = {}
    dataset_annotation_names = []
    dataset_names = []
    for dataset in datasets_yaml:
        dataset_names.append(dataset)
        dataset_annotation_names.append(datasets_yaml[dataset][yaml_keys.kDAnnotation])
#         dataset_annotation_names_dict[datasets_yaml[dataset][yaml_keys.kDAnnotation]] = 1
    datasets_db = db_utils.get_table_names(c)
    # check db against datasets.yml
    for dataset in datasets_db:
        if re.sub('_r$', '', dataset) not in dataset_annotation_names:
            warning = 'warning: ' + str(dataset) + ' in db but not datasets.yml (or commented out)'
            warnings.append(warning)
            print warning
    # check datasets.yml against db
    for idx,dataset in enumerate(dataset_annotation_names):
        dataset_yaml_name = dataset_names[idx]
        if(dataset not in datasets_db and dataset+'_r' not in datasets_db):
            if(datasets_yaml[dataset_yaml_name][yaml_keys.kDImportIfMissing]):
                error = 'error: ' + str(dataset) + ' in datasets.yml but not database. Run stmp.py --update_db to load it into the database, or comment it out (using # before each line) or change Import_If_Missing to False in datasets.yml to ignore it.'
                errors.append(error)
                print error
            elif(not datasets_yaml[dataset_yaml_name][yaml_keys.kDImportIfMissing]):
                warning = 'warning: ' + str(dataset) + ' in datasets.yml but not database, and is currently not being imported (Import_If_Missing = False). If desired, change Import_If_Missing to True for this dataset in datasets.yml and run stmp.py --update_db to load it into the database. Otherwise this dataset will not be used for annotations.'
                warnings.append(warning)
                print warning
        
    # 2-2) check dataset columns in yaml vs db (both directions)
    # WARNING: we do not currently check the delimiter for multiple matches, if stored in the db (specified in the defaults section of datasets.yml)
    for dataset in datasets_yaml:
        dataset_annotated_name = datasets_yaml[dataset][yaml_keys.kDAnnotation]
        dataset_db_name = ''
        if(dataset_annotated_name in datasets_db):
            dataset_db_name = dataset_annotated_name
        elif(dataset_annotated_name+'_r' in datasets_db):
            dataset_db_name=dataset_annotated_name+'_r'
        else: # dataset not found in db (which means we already gave an error/warning above)
            continue
        if(db_utils.is_region_dataset(dataset_db_name) != yaml_utils.is_region_dataset(dataset, yaml_commands)):
            error = 'error: dataset ' + str(dataset_db_name) + ' in db does not match with ' + str(dataset) + ' in datasets.yml. One is specified as region and the other is specified as point. Either change datasets.yml or remove and reimport the dataset into the database.'
            errors.append(error)
            print error
        # 2-2a) check column names in yaml vs db
        yaml_cols = datasets_yaml[dataset][yaml_keys.kDColumnHeaders]
        db_cols = db_utils.getColumnNamesOfTable(c, dataset_db_name)
        for col in yaml_cols:
            if(col != '' and col not in db_cols 
               and col not in stmp_annotation_util.START_COLUMN_HEADERS and col not in stmp_annotation_util.STOP_COLUMN_HEADERS and col not in stmp_annotation_util.CHR_HEADERS):
                # TODO check for whether col is start, stop, or chr col and then check whether these cols are in the imported table with standard names (for true db integrity check)
                error = 'error: column ' + str(col) + ' of dataset ' + str(dataset) + ' found in datasets.yml but not database'
                errors.append(error)
                print error
        # 2-2b) check column names in db vs yaml
        for col in db_cols:
            if(col not in yaml_cols 
               and col.lower() not in stmp_annotation_util.START_COLUMN_HEADERS and col.lower() not in stmp_annotation_util.STOP_COLUMN_HEADERS and col.lower() not in stmp_annotation_util.CHR_HEADERS 
               and col.lower() not in vcfHeaders.kVCFColTypes 
               and col != vcfHeaders.kClinvarStarHeader):
                error = 'error: column ' + str(col) + ' of dataset ' + str(dataset) + ' found in database but not datasets.yml (and is not a standard VCF column, clinvar star rating, chr, start, or stop column)'
                errors.append(error)
                print error
    
    # (End of checks) summarize with all the errors, warnings, etc.
    err_warning_summary = 'Warnings:\n' + '\n'.join(warnings) + '\n\n' + 'Errors:\n' + "\n".join(errors) + '\n\n' + str(len(warnings)) + ' warnings and ' + str(len(errors)) + ' errors in preflight checks.'
    if(len(errors) > 0):
        if(raise_exception):
            raise ValueError(err_warning_summary + ' Please correct the errors and then re-run STMP.')
        #else
        return False
    #else
    print "\n" + err_warning_summary
    print '+(OK) passed pre-flight checks\n'
    return True
    
