#!/usr/bin/env python

'''
Author: Prag Batra prag@stanford.edu

Purpose:
    General helper functions for STMP (e.g. resolving relative vs absolute paths, etc.)
    

Explanation:

    

Example:

    

'''

import os

def list2dict(list):
    dict = {}
    for elt in list:
        dict[elt] = 1
    
    return dict

# assumes headlist is first line of file
def get_headlist_tsv(myfile):
    f = open(myfile, 'r')
    lineContents = f.readline().rstrip("\n").split("\t")
    return lineContents

# if needed, converts relative path (relative to the current working directory) to an absolute path
def root_or_cwd(mydir):
    # If the user has specified a file using a partial filepath relative to the current working directory, complete that path so that 
    # the file may be located.  If the filepath is absolute (i.e. starting from the root directory) then leave it alone.
    
    if mydir[0] != '/':
        return os.path.join(os.getcwd(), mydir) # merge the current working directory to the provided partial path
    else:
        return(mydir) # no change

# gets name of file/directory (including extension)
def get_file_or_dir_name(path):
    return os.path.split(path)[-1]

