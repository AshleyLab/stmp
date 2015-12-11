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

The script does so in a streaming fashion from stdin.
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
	if(line.startswith('#')): # header
		print (line.rstrip("\n"))
		continue
	#else
	s = line.rstrip('\n').split("\t")
	if firstTime:
# 		prevLine = ['']*len(s)
		prevLine = s
		firstTime = False
		continue
	#else
	# 3 cases: 
	# 1) duplicate annotation (same chr, start, stop, annotation)
	# 2) additional annotation (same chr, start, stop (aka same variant), different annotation)
	# 3) different annotation (different chr, start, and/or stop, aka different variant)
	if(not firstTime):
		# case 1: duplicate annotation
		if(s[0] == prevLine[0] and s[1] == prevLine[1] and s[2] == prevLine[2] and s[3] == prevLine[3]): # chr, start, stop, annotation match
			continue # do nothing -- we already have this annotation in the annotations list
		# case 2: additional annotation
		elif(s[0] == prevLine[0] and s[1] == prevLine[1] and s[2] == prevLine[2] and s[3] != prevLine[3]): # chr, start, stop match, but not annotation (additional annotation)
			# add to list of annotations for this variant
			annotations.append(s[-1])
		# case 3: different annotation (annotation for a different variant)
		# trigger output + reset annotation list
		else:
			print ('\t'.join(prevLine) + '\t' + bed_multimatch_internal_delimiter.join(annotations)) # for now, print out entire VCF line + annottion (TODO print out just annotation in future)
			annotations = []
			annotations.append(s[-1])
			
	prevLine = s

# print out any remaining annotations
if(len(annotations) > 0):
	print ('\t'.join(prevLine) + '\t' + bed_multimatch_internal_delimiter.join(annotations)) # for now, print out entire VCF line + annotation (TODO print out just annotation in future)
