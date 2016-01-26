#!/usr/bin/env python

'''
Author: Prag Batra prag@stanford.edu

Purpose:
    General helper functions for STMP (e.g. resolving relative vs absolute paths, etc.)
    

Explanation:

    

Example:

    

'''

import os
import datetime
import subprocess
import gzip
import mmap


# improved open function that handles .gz files properly
def open_compressed_or_regular(f, options):
    if(f.endswith('.gz')):
        return gzip.open(f, options)
    #else
    return open(f, options)

# gets number of lines in a file
def get_num_lines(filePath, ignore_vcf_info_lines=True):
    if(ignore_vcf_info_lines):
        cmd = 'grep -v "##" {file}|wc -l'.format(file=filePath)
    else:
        cmd = 'wc -l {file}'.format(file = filePath)
#     cmd = 'wc -l {file}'.format(file=filePath)
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    stdout,stderr = proc.communicate()
    returncode = proc.returncode
    if(returncode != 0):
        raise ValueError('Failed to get number of lines in file. Cmd: ' + cmd + '\nError: ' + str(stderr) + '\nOutput: ' + str(stdout) + '\nReturn Code: ' + str(returncode))
    #else
    numLines = int(stdout.split(' ')[0])
    return numLines


# gets absolute path to directory where the code (i.e. this script file) is
def get_code_dir_abs():
    return os.path.dirname(os.path.realpath(__file__))

# prints the given string every [interval] s or longer (depending on how often this method is called). Call it each time the progress (e.g. line number) changes.
# NOTE: prevTime and interval must be in SECONDS
def print_progress_timed(progressStr, prevTime, interval=5):
    currtime = datetime.datetime.now().second
    if(abs(currtime - prevTime) > interval):
        print str(datetime.datetime.now()) + ': ' + progressStr
        return currtime
    #else
    return prevTime

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

# if needed, converts relative path to absolute path (with respect to code directory)
def root_or_code_dir(mydir):
    if mydir[0] != '/':
        return os.path.join(get_code_dir_abs(), mydir)
    else:
        return mydir

# gets the absolute path to the parent dir of the given path
def get_parent_dir(path):
    return os.path.abspath(os.path.join(path, os.pardir))

# gets name of file/directory (including extension)
def get_file_or_dir_name(path):
    return os.path.split(path)[-1]


# CURRENTLY UNUSED
# helper function: uses mmap to delete a given set of positions from a file
def deleteFromMmap(f, mm, start, end):
    length = end - start
    size = len(mm)
    newsize = size - length

    mm.move(start,end,size-end)
    mm.flush()
    mm.close()
    f.truncate(newsize)
    mm = mmap.mmap(f.fileno(),0)
    return mm

