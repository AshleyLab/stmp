#!/usr/bin/env python

'''
Author: Prag Batra prag@stanford.edu

Purpose:
    Helper methods for parsing YAML configuration files (e.g. for STMP)
    

Explanation:

    

Example:

    

'''
import yaml
import yaml_keys
import os

def convertColumns(cols, yaml_commands):
    actualCols = []
    for col in cols:
        if(isinstance(col, list)):
            convertedColList = convertColumns(col, yaml_commands) # recursion
            actualCols.append(convertedColList)
            continue
        if('.' not in col): # nothing to do here
            actualCols.append(col)
            continue
        #else
        colComponents = col.split('.')
        dataset_yaml = yaml_commands[colComponents[0]]
        dataset_annotation_name = dataset_yaml[yaml_keys.kDAnnotation]
        if(dataset_yaml[yaml_keys.kDCategory] == yaml_keys.kDCategoryTypeRegion):
            dataset_annotation_name += '_r'
        actualCols.append(dataset_annotation_name+'_'+colComponents[1])

    return actualCols

#converts tiering column references to actual headers
# def convertTieringColumns(yaml_commands):
#     return convertColumns(yaml_commands[yaml_keys.kModules][yaml_keys.kTiering], yaml_commands)
#     
#     #debug
# #     print 'yaml_commands[{modules} {tiering} {alleleFreq}'.format(modules=yaml_keys.kModules, tiering=yaml_keys.kTiering, alleleFreq=yaml_keys.kTAlleleFreqCols)
#     
#     allele_freq_cols = yaml_commands[yaml_keys.kModules][yaml_keys.kTiering][yaml_keys.kTAlleleFreqCols]
#     actualCols = []
#     for freqCol in allele_freq_cols:
#         freqColComponents = freqCol.split('.')
#         actualCols.append(yaml_commands[freqColComponents[0]][yaml_keys.kDAnnotation]+'_'+freqColComponents[1])
#     
# #     yaml_commands[yaml_keys.kTiering][yaml_keys.kTAlleleFreqCols] = 
#     return actualCols

# loads YAML file and returns content as nested python dictionaries/arrays
def parse_yaml(loc):
    # Parse a YAML file which will instruct the annotation steps to be completed.
    with open(loc, "r") as stream:
        yaml_commands = yaml.safe_load(stream)
        return yaml_commands

# combines yaml input files, parses, and outputs single python structure containing commands (nested dictionaries, etc.)
def parse_yaml_input_files(config_file, modules_file):
    config_cmds = parse_yaml(config_file)
    module_cmds = parse_yaml(modules_file)
    config_cmds[yaml_keys.kModules] = module_cmds
    return config_cmds

# gets absolute path wrt location of this file (which should be in the same directory as stmp.py)
def get_abs_path(yaml_path):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    if(not yaml_path.startswith('/')): # not an absolute path
        yaml_path = os.path.join(script_dir, yaml_path)
    return yaml_path

# splits a given set of threshold values for a column using the separator specified in the YAML tiering config
def split_multiple_col_thresholds(col_threshold_str, yaml_commands):
    col_threshold_separator = yaml_commands[yaml_keys.kModules][yaml_keys.kTiering][yaml_keys.kTColMultipleThresholdSeparator]
    return col_threshold_str.split(col_threshold_separator)
