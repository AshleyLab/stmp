#!/usr/bin/env python

'''
Author: Prag Batra prag@stanford.edu

Purpose:
    Expands intersectBed output files so they have the appropriate number of columns (for datasets with multiple annotation columns). Uses the delimiters specified in modules.yml and datasets.yml.
    

Explanation:

    

Example:

    

'''

import sys
from sys import argv
import argparse
import yaml_utils
import yaml_keys
import general_utils
import datetime


def expandBED(bedPath, bed_multimatch_internal_delimiter, bed_delimiter, dataset_multimatch_delimiter, outSuffix):
    print 'Expanding file ' + str(bedPath)
    f = open(bedPath, 'r')
    fileHeader = f.readline().rstrip("\n").split("\t")
    f.seek(0) # reset file position so we can reread the header below
    out_file_path = bedPath+outSuffix
    out_file = open(out_file_path, 'w')
    prevTime = datetime.datetime.now().second
    for currLine,fl in enumerate(f):
        prevTime = general_utils.print_progress_timed('expanding line ' + str(currLine), prevTime)
        lineContents = []
        fl = fl.rstrip("\n")
        if(fl == '.' or len(fl) == 0):
            # extend lineContents by appropriate # of blank columns
            lineContents.extend(['']*len(fileHeader))
        elif(bed_multimatch_internal_delimiter in fl):
            # split annotations and then recombine using specified default delimiter
            annotations = fl.split(bed_multimatch_internal_delimiter)
#             finalColsDict = [{}]*len(fileHeader) # array of dicts, one for each col
            finalCols = ['']*len(fileHeader) # ensures the appropriate number of cols for future hits (which will be appended to the appropriate empty string)
            annotationsDict = {}
            for annotation in annotations:
                if(annotation in annotationsDict):
                    continue
                #else
                annotationsDict[annotation] = 1
                annotationComponents = annotation.split(bed_delimiter)
                if(len(annotationComponents) > 0):
                    for idx,component in enumerate(annotationComponents):
                        # WARNING: we are not checking alignment of annotations across columns here, so if there are some columns within the dataset that have multiple annotations while others are completely blank, a given annotation in one column may not correspond with the annotation at the same position in the other column (e.g. annotation#2 in column 1 may not correspond to annotation#2 in column 2).
                        try:
                            if(finalCols[idx] != ''):
                                finalCols[idx] += dataset_multimatch_delimiter
                            finalCols[idx] += component
                        except IndexError:
                            print 'error: index ' + str(idx) + ' beyond size of finalCols (' + str(len(finalCols)) + ') on file ' + str(bedPath) + ' line ' + str(currLine)
                            raise 
            
#             # run thru final cols again and remove duplicate entries
#             for col in finalCols:
#                 
            lineContents.extend(finalCols)
        else:
            flc = fl.split(bed_delimiter)
            lineContents.extend(flc)
        out_file.write("\t".join(lineContents) + "\n")
    return out_file_path


##### MAIN CODE ####
if __name__ == '__main__':
    # parse commandline args.
    parser = argparse.ArgumentParser(description = 'IntersectBed Expander')
    
    parser.add_argument('--tab_file', help='.tab file to be expanded')
    parser.add_argument('--out_suffix', default='.tsv', help='suffix for output file (stored in same location as input file), e.g. ".tsv"')
    parser.add_argument('--config_modules', help='path to modules configuration file (YAML)')
    parser.add_argument('--config_datasets', help='path to datasets configuration file (YAML)')
    
    args = parser.parse_args()
    
    config_yaml = yaml_utils.parse_yaml_input_files(args.config_datasets, args.config_modules)
    modules_yaml = config_yaml[yaml_keys.kModules]
#     modules_yaml = yaml_utils.parse_yaml(args.config_modules)
#     datasets_yaml = yaml_utils.parse_yaml(args.config_datasets)
    
    bed_delimiter = modules_yaml[yaml_keys.kAnnotation][yaml_keys.kABedInternalDelimiter]
    bed_multimatch_internal_delimiter = modules_yaml[yaml_keys.kAnnotation][yaml_keys.kABedMultimatchInternalDelimiter]
    # for now, just using the default delimiter. Later, check if there's a specific delimiter for a given dataset and use that instead.
    dataset_multimatch_delimiter = yaml_utils.get_dataset_defaults(config_yaml)[yaml_keys.kDMultimatchDelimiter]
    
    expandBED(args.tab_file, bed_multimatch_internal_delimiter, bed_delimiter, dataset_multimatch_delimiter, args.out_suffix)

