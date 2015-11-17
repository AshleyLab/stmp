#!/usr/bin/env python
'''
Author: Eli Moss elimoss@stanford.edu, Prag Batra prag@stanford.edu


Purpose: 

	Modify the output of bedtools' intersectbed utility in order to contain only one instance of each genomic locus.  

Explanation: 

	The left outer join functionality of intersectbed will output one line per match between its two inputs, whereas we really want only one line 
	per genomic locus.  The purpose of this script then is to take the output of intersectbed and, when a locus appears more than once
	with multiple annotations, produce one line for that locus with the annotations joined with commas.

Example:

	1	241320	rs79974410	A	G	3.22	PASS	DP	70	13_Heterochrom/lo
	1	241320	rs79974410	A	G	3.22	PASS	DP	70	11_Weak_Txn

	becomes

	1	241320	rs79974410	A	G	3.22	PASS	DP	70	13_Heterochrom/lo,11_Weak_Txn

The script does in a streaming fashion from stdin.
'''

from __future__ import print_function
import sys
import argparse
import yaml_keys
import yaml_utils

parser = argparse.ArgumentParser()
parser.add_argument('--modules', help='modules yaml config file (to get BED delimiter info)')
args = parser.parse_args()
modules = yaml_utils.parse_yaml(args.modules)
bed_multimatch_internal_delimiter = modules[yaml_keys.kAnnotation][yaml_keys.kABedMultimatchInternalDelimiter]

annotations = []

firstTime = True

for line in sys.stdin:
	s = line.rstrip('\n').split("\t")
	if firstTime:
		prevline = [''] * len(s)
# 	print ('number of entries: ' + str(len(s)), file=sys.stderr) #debug
# 	print ('s[1] ' + s[1] + ' pl[1] ' + prevline[1] + ' s2 ' + s[2] + ' pl2 ' + prevline[2] 
# 	+ ' s3 ' + s[3] + ' pl3 ' + prevline[3], file=sys.stderr) #debug
	if not firstTime and (line.startswith('#') or s[0] != prevline[0] or s[1] != prevline[1] or s[2] != prevline[2] or s[3] != prevline[3] or s[4] != prevline[4]): # line doesn't match previous, or it's a header.  In either case, trigger output.
		print('\t'.join(prevline[0:12] + [bed_multimatch_internal_delimiter.join(annotations)]))
		annotations = []
	else: # it's a duplicate!  Start/continue annotation collection.
		annotations.append(s[-1])
		firstTime = False

	prevline = s

# print out any remaining annotations
if len(annotations) > 0:
	print('\t'.join(prevline[0:12] + [bed_multimatch_internal_delimiter.join(annotations)]))

