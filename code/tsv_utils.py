#!/usr/bin/env python

'''
Author: Prag Batra prag@stanford.edu

Purpose:
    Utilities (helper functions) for analyzing and manipulating annotated TSV files (e.g. those output by STMP annotation module).
    

Explanation:

    

Example:

    

'''

# Gets TSV header as a dictionary (key = header of a given column, value = index of that column)
# NOTE: file handle must be at beginning of file!! (assumes 1st row is header with no #)
# NOTE: does not reset position in file! Use seek() if you want to go back to beginning of file.
def getTSVHeaderCols(fileHandle):
    headerStr = fileHandle.readline()
    header = headerStr.split("\t")
    headerDict = {}
    for idx,h in enumerate(header):
        headerDict[h] = idx
    return headerDict

