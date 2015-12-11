#!/usr/bin/env python

'''
Author: Prag Batra prag@stanford.edu

Purpose:
    Performs basic checks on STMP annotated output file (e.g. no columns that are completely empty, etc.)
    

Explanation:

    

Example:

    

'''

import sys
from sys import argv
import argparse
import copy
import subprocess

import general_utils

# checks files in parallel to make sure they all have the same line count. By default, raises
# an exception if this is not the case. If raise_exception is set to false, prints a warning instead. If lineCount is left at the default (-1), uses the number of lines in the first specified file for comparison against all other files.
def check_file_numlines(files, lineCount=-1, ignore_vcf_info_lines=True, raise_exception=True):
    bedPaths = files
    lineCountBedFile = ''
    retVal = True
    for bedFile in bedPaths:
        if(ignore_vcf_info_lines):
            cmd = 'grep -v "##" {file}|wc -l'.format(file=bedFile)
        else:
            cmd = 'wc -l {bedFile}'.format(bedFile = bedFile)
        print 'cmd to check # of lines in file: ' + cmd
        processes = []
        processes.append(subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE))
    for idx,process in enumerate(processes):
        out,err = process.communicate() # blocks until process terminates
        returnvalue = process.returncode
        if(returnvalue != 0):
            raise ValueError('Failed to get number of lines for file. Cmd: ' + str(cmd) + '\nError: ' + str(err) + '\nOutput: ' + str(out) + '\nReturn code: ' + str(returnvalue))
        #else
        outarr = out.split(' ')
        #debug
#         print 'cmd output as array: ' + str(outarr)
        numLines = int(outarr[-1])
        if(lineCount == -1):
            lineCount = numLines
            lineCountBedFile = bedPaths[idx]
        elif(numLines != lineCount):
            if(lineCountBedFile != ''):
                errorText = str(numLines) + ' lines in ' + str(bedPaths[idx]) + ' different from expected ' + str(lineCount) + ' lines in ' + str(lineCountBedFile) + '\ncmd output as array: ' + str(outarr)
            else:
                errorText = str(numLines) + ' lines in ' + str(bedPaths[idx]) + ' different from expected ' + str(lineCount) + ' lines' + '\ncmd output as array: ' + str(outarr)
            
            if(raise_exception):
                raise ValueError(errorText)
            else:
                print 'Warning: ' + errorText
                retVal = False
    
    if(not retVal):
        return False
    #else
    print '+(OK) Files all have same line count (' + str(lineCount) + ')'
    return True
    
    
# currently checks annotated output for any columns that are completely empty (no value for any row in the output file)
def check_annotated_output(out_file, suppress_print = False): 
    f = open(out_file, 'r')
    headlist = f.readline().rstrip("\n").split("\t")
    
    filledColsDict = {}
    for colHeader in headlist:
        filledColsDict[colHeader] = 0
    
    for line in f:
        lineContents = line.rstrip("\n").split("\t")
        for idx,col in enumerate(lineContents):
            if col != '' and col != ' ':
                colHeader = headlist[idx]
                filledColsDict[colHeader] = 1 # not counting, just simple binary check
    
    emptyCols = []
    
    for colHeader in filledColsDict:
        if filledColsDict[colHeader] == 0:
            emptyCols.append(colHeader)
    
    if(not suppress_print):
        print 'Warning: the following {num} of {totalNum} annotated output columns have no value for any variant in this file: '.format(num=str(len(emptyCols)), totalNum=str(len(headlist)))
        for emptyCol in emptyCols:
            print emptyCol
        
    return emptyCols

# can also be run from command line
if __name__ == '__main__':
    # parse commandline args.
    parser = argparse.ArgumentParser(description = 'STMP Annotation Checker')
    parser.add_argument('--input_files', help='comma-separted list of input files to check')
    
    args = parser.parse_args()
    
    #check header for each file to see if same set of cols is present in each file
    global_cols = {}
    file_header_dicts = {}
    
    input_files = args.input_files.split(",")
    for input_file in input_files:
        cols = general_utils.get_headlist_tsv(input_file)
        file_header_dicts[input_file] = general_utils.list2dict(cols)
        for col in cols:
            global_cols[col] = 1
    
    for file in file_header_dicts:
        headerDict = file_header_dicts[file]
        for col in global_cols:
            if(col not in headerDict):
                print 'warning: col ' + str(col) + ' not in ' + str(file)
    
    overall_missing_cols = copy.deepcopy(global_cols)
    
    for input_file in input_files:
        missing_cols_dict = general_utils.list2dict(check_annotated_output(input_file, suppress_print=True))
        
        new_missing_cols = copy.deepcopy(overall_missing_cols)
        for col in overall_missing_cols:
            if(col not in missing_cols_dict):
                del new_missing_cols[col]
        overall_missing_cols = new_missing_cols
    
    print str(len(overall_missing_cols)) + ' out of {numcols} cols with no value across ANY of the {num} input files: '.format(num=str(len(input_files)), numcols=len(global_cols))
    for col in overall_missing_cols:
        print str(col)

