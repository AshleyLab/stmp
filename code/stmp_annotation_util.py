#!/usr/bin/env python

'''
Author: Eli Moss elimoss@stanford.edu, Prag Batra prag@stanford.edu

Purpose:
	
	All the heavy lifting for variant annotation happens here.

Explanation:

	Functionality is so diverse that descriptions are best left to inline comments.

Example:

	For usage, see stmp.py where these methods are all called.

'''

import yaml_keys
import sys
import linecache
import sqlite3
import os
from subprocess import *
import subprocess
import stmp_consts
import yaml
import mmap
import datetime
import vcfHeaders
import vcfUtils
import tsv_utils
import gzip
import yaml_utils
import general_utils
import multiprocessing
from multiprocessing import Pool
import expand_intersectedBeds
import stmp_annotation_checker
import db_utils

# UNUSED - Locations of important things for snpeff (ideally should be in user's path). For now we're specifying them in modules.yml instead.
# SNPEFF = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, 'third_party', 'snpeff', 'snpEff', 'snpEff.jar'))
# SNPEFF_MEMORY_ALLOCATION = '6g'

# still used but should be switched to modules.yml soon.
ANNOVAR_PATH = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, 'third_party', 'annovar')) # we need this both for annovar executable (annotate_variation.pl, table_annovar.pl) and to locate humandb folder for functional annotations

# still used but should be switched to modules.yml soon.
START_COLUMN_HEADERS = ['start', 'begin', 'pos']
STOP_COLUMN_HEADERS = ['stop', 'end']
CHR_HEADERS = ['chr', 'chrom', 'chromosome']

# convenience method for getting annovar path from YAML (currently unused)
def get_annovar_path(yaml_commands):
	return os.path.join(general_utils.get_code_dir_abs(), yaml_commands[yaml_keys.kModules][yaml_keys.kAnnotation][yaml_keys.kAAnnovarPath])


def open_compressed_or_regular(f, options):
	return general_utils.open_compressed_or_regular(f, options)

def downloadDB(db_yaml, out_path, out_dir, download_method, ucsc_basePath = None):
	if(download_method == 'ANNOVAR'):
		# TODO use annovar to download
		print 'Error ANNOVAR download method not yet implemented for dataset ' + db_yaml[yaml_keys.kDAnnotation]
		None
	elif(download_method == 'URL'):
		if not os.path.exists(out_dir):
			os.makedirs(out_dir)
		cmd = "cd {out_dir}; wget --retr-symlinks --no-check-certificate --timestamping {url}".format(out_dir=out_dir, outfile=out_path, url=db_yaml[yaml_keys.kDSource]) # For now, just download the file with the same name as on the server (we'll use the name given in the YAML when importing it into our database). "--timestamping" flag won't re-retrieve files unless newer than local file.
		print 'command to download dataset: ' + str(cmd)
		subprocess.Popen(cmd, shell=True).wait()
		#record download date in separate metadata file for versioning purposes
		metadata_path = out_path+'.info'
		writef(metadata_path, 'Downloaded '+str(datetime.datetime.now())+"\n")  #TODO maybe capture download time at beginning insead of end of download
		
		return out_path
		# TODO maybe include UCSC stuff as part of URL download (instead of having UCSC base path, etc.)
		None
# 	elif(download_method == 'UCSC'):
# 		# TODO use UCSC base path to download dataset
		

# downloads datasets from annovar (via annovar or UCSC website)
def downloadDBs(yaml_commands, out_dir, logfile = None):
	yaml_datasets = yaml_utils.get_datasets(yaml_commands)
	yaml_dataset_defaults = yaml_utils.get_dataset_defaults(yaml_commands)
	out_log = None
	if logfile != None:
		out_log = open(logfile, 'w')
	
	build = yaml_dataset_defaults[yaml_keys.kDBuild]
	
	for dataset in yaml_datasets:
# 		if(dataset == 'default'):
# 			build = yaml_commands[dataset]['Build']
# 			continue
# 		if(dataset==yaml_keys.kModules): # skip module spec as that isn't a database
# 			continue
		
		if(yaml_datasets[dataset] != None and yaml_keys.kDAnnotation in yaml_datasets[dataset]):
			dataset_name = yaml_datasets[dataset][yaml_keys.kDAnnotation]
			print 'Downloading dataset for annotation: ' + dataset_name
			# attempt to download from annovar site
			cmd = 'perl {annotate_variation} -buildver {build} -downdb -webfrom annovar {dbname} {out_dir}'.format(build=build, annotate_variation=os.path.join(ANNOVAR_PATH, 'annotate_variation.pl'), dbname=dataset_name, out_dir=out_dir)
			print 'annovar db download command: ' + str(cmd)
			retVal = subprocess.Popen(cmd, shell=True).wait()
			# TODO grab name/path of annovar downloaded file
			if(retVal != 0):
				# attempt to download from UCSC site
				cmd2 = 'perl {annotate_variation} -buildver hg19 -downdb {dbname} {out_dir}'.format(annotate_variation=os.path.join(ANNOVAR_PATH, 'annotate_variation.pl'), dbname=dataset_name, out_dir=out_dir)
				print 'Retrying with Annovar alt db download command: ' + str(cmd2)
				retVal = subprocess.Popen(cmd2, shell=True).wait()
				# TODO grab name/path of downloaded file
				if(retVal != 0):
					print 'unable to download ' + dataset_name
					if out_log != None:
						out_log.write('unable to download ' + dataset_name + "\n")
		
	print 'Finished downloading datasets to ' + out_dir


# checks which DBs were downloaded
def checkDBs(yaml_commands, out_dir):
	print 'checking downloaded datasets at ' + out_dir + ' against YAML'
	files = os.listdir(out_dir)
	filesHash = {}
	for file in files:
		filesHash[os.path.splitext(file)[0].lower()] = 1
	
	missingDatasets = 0
	yaml_datasets = yaml_utils.get_datasets(yaml_commands)
	for dataset in yaml_datasets:
		if('hg19_' + dataset.lower() not in filesHash):
			print 'missing dataset ' + dataset
			missingDatasets += 1
	
	print 'Total datasets in YAML: ' + str(len(yaml_datasets.keys()))
	print 'Missing datasets: ' + str(missingDatasets)


def getSampleTableName(vcf_file_loc):
	return 'sample_' + getSampleName(vcf_file_loc)

def getSampleName(vcf_file_loc):
	return vcf_file_loc[0:len(vcf_file_loc) - 4].replace('.', '_').replace('-', '_').split('/')[-1]

def getFilename(file_path):
	return os.path.split(file_path)[-1]

def stripChrPrefix(vcf_filepath, out_dir, skip_if_exists=False):
	print 'Stripping chr prefix from VCF CHROM column (if present)'
	vcf_name = getSampleName(vcf_filepath)
	outfilepath = os.path.join(out_dir, vcf_name+'_strippedChr' + '.vcf.gz')
	if(skip_if_exists and os.path.isfile(outfilepath)):
		print 'Reusing existing stripped chr VCF file at ' + outfilepath
		return outfilepath
	#else
	#strip chr prefix
	vcf = open_compressed_or_regular(vcf_filepath, 'r')
	out = open_compressed_or_regular(outfilepath, 'w')
	
	for line in vcf:
		if(line.startswith('#')):
			out.write(line)
			continue
		#else
		lineContents = line.rstrip("\n").split("\t")
		chrom = lineContents[0]
		if(chrom.startswith('chr')):
			chrom = chrom[3:] # strip 'chr' prefix
		lineContents[0] = chrom
		newline = "\t".join(lineContents)+"\n"
		out.write(newline)
	
	return outfilepath


# splits multiallelic site in VCF into separate rows using bcftools
def splitMultiallelic(vcf_filepath, out_dir, skip_if_exists=False):
	print 'Converting VCF so there is only 1 allele per line'
	vcf_name = getSampleName(vcf_filepath)
	outfilepath = os.path.join(out_dir, vcf_name + '_noMultiallelic' + '.vcf.gz')
	if(skip_if_exists and os.path.isfile(outfilepath)):
		print 'Reusing existing converted VCF file at ' + outfilepath
		return outfilepath
	#else
	cmd = "bcftools norm -m - '{vcf_filepath}' -O z -o '{out_dir}'".format(vcf_filepath=vcf_filepath, out_dir=outfilepath)
	print 'bcftools split multiallelic command: ' + str(cmd)
	
	subprocess.Popen(cmd, shell=True).wait()
	return outfilepath


def stripVCFComments(vcf_filepath, out_dir):
	print 'Removing extra fields (##) from VCF'
	vcf_name = getSampleName(vcf_filepath)
	outfilepath = os.path.join(out_dir, vcf_name + '_noComments' + '.vcf')
	cmd = 'grep -v "^##" {vcf_filepath} > {vcf_outfile}'.format(vcf_filepath = vcf_filepath, vcf_outfile = outfilepath)
	print 'command to strip vcf comments: ' + str(cmd)
	subprocess.Popen(cmd, shell=True).wait()
	return outfilepath


def annovar_annotate_functional(vcf_file_loc, output_dir, buildver='hg19', humandb_loc=os.path.join(ANNOVAR_PATH, 'humandb'), genedefs=['refGene', 'knownGene', 'wgEncodeGencodeBasicV19']):
	sample_name = getSampleName(vcf_file_loc)
	print('Running annovar on ' + sample_name)
	
	outpathprefix = os.path.join(output_dir, 'annovar')
	
	# annotate
	cmd = '{table_annovar} {vcf_file_loc} {humandb_loc} -buildver {buildver} -out {outprefix} -remove -protocol {protocols} -operation {operations} -vcfinput'.format(vcf_file_loc=vcf_file_loc, table_annovar = os.path.join(ANNOVAR_PATH, 'table_annovar.pl'), humandb_loc=humandb_loc, buildver=buildver, outprefix=outpathprefix, protocols=",".join(genedefs), operations=','.join("g" for i in range(0, len(genedefs))))
	
	print 'annovar annotation cmd: ' + str(cmd)
	
	outfile = os.path.join(output_dir, 'annovar' + '.' + buildver + '_multianno.txt')
	
	return [subprocess.Popen(cmd, shell=True), outfile]


# helper function to convert a VCF consisting of per-sample allele info to a VCF with just allele freqs using bcftools
# currently waits for completion, but could be easily made async (immediately returning what the output file name would be after completion)
def vcf2allelefreq(vcf_file_loc, output_dir, overwriteIfExists=False):
	sample_name = getSampleName(vcf_file_loc)
	outfilename = sample_name + '_alleleFreqs.vcf'
	outfilepath = os.path.join(output_dir, outfilename)
	if((not overwriteIfExists) and os.path.isfile(outfilepath)):
		print 'skipped conversion of ' + sample_name + ' to allele freq VCF because already converted'
	else:
		cmd = "bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC\t%AN\n' " + vcf_file_loc + " | awk '{OFS=\"\t\";print $0,$6/$7}' > " + outfilepath
		print 'vcf to allele freq bcftools cmd: ' + str(cmd)
		subprocess.Popen(cmd, shell=True).wait()
	
	return outfilename

def snpeff(vcf_file_loc, output_dir, yaml_commands):
	# In a subprocess, start a SnpEff run on the sample VCF, and output it to the user-designated folder.
	# By returning the process pointer, a step where its completion is awaited can be added to the top-level
	# master script.  That way we don't need to hold things waiting for completion here.
	snpeff_memory_allocation = yaml_commands[yaml_keys.kModules][yaml_keys.kAnnotation][yaml_keys.kASnpeffMemory]
	snpeff_jar_path = os.path.join(general_utils.get_code_dir_abs(), yaml_commands[yaml_keys.kModules][yaml_keys.kAnnotation][yaml_keys.kASnpeffPath], 'snpEff.jar')
	sample_name = getSampleName(vcf_file_loc)
	print("Running snpEff on " + sample_name)
	cmd = "java -Xmx{memalloc} -jar {snpeff_jar} eff hg19 {vcf} -stats {out_loc} -csvStats".format(snpeff_jar = snpeff_jar_path, vcf = vcf_file_loc, memalloc=snpeff_memory_allocation, out_loc=output_dir)
	#debug
	print 'snpeff annot cmd: ' + cmd
	annotated_vcf = os.path.join(output_dir, "{sample_name}.snpeff.vcf".format(sample_name = sample_name))
	p = subprocess.Popen(cmd, stdout = open_compressed_or_regular(annotated_vcf, 'w'), shell = True)
	
	return [p, sample_name, annotated_vcf]


# converts snpeff vcf file to tsv by removing lines with '#'
# snpeff_vcf_outfile = directory of this file
def snpeff2tsv(sample_name, snpeff_vcf_outfile, output_dir):
	# strip comments from VCF to make it a tsv
	cmd = "grep -v '^##' {vcf}".format(vcf = snpeff_vcf_outfile)
	print 'snpeff vcf to tsv conversion command (remove \'##\' lines): ' + str(cmd)
	tsv_out_filePath = os.path.join(output_dir, "{sample_name}.snpeff.tsv".format(sample_name = sample_name))
	out_file = open(tsv_out_filePath, 'w')
	subprocess.Popen(cmd, stdout = out_file, shell=True).wait()
	
	return tsv_out_filePath

#writes the given text to the file, overwriting anything else in the file
# NOTE: automatically writes "\n" at end of text (similar to print command)
def writef(file_path, text):
	f = open_compressed_or_regular(file_path, 'w')
	f.write(text + "\n")
	f.close()


def getVCFHeader(vcf_file_loc):
	vcf = open_compressed_or_regular(vcf_file_loc, 'r')
	for line in vcf:
		if(line[0] == '#' and line[1] != '#'):
			vcf.close()
			return line[1:].rstrip("\n")
	
	# else
	print 'Error no header found in vcf ' + vcf_file_loc
	vcf.close()
	return ''

#returns filename without extension
def stripExtension(filename):
	return filename.split('.')[0]

# converts VCF to TSV, extracting each specified tag from the INFO column and listing it as a seperate column in the output TSV
def vcf2tsv(vcf_path, tags, output_dir, output_extension='.tsv'):
	vcf_filename = getFilename(vcf_path)
	vcf_filename_noExtension = stripExtension(vcf_filename)
	outpath = os.path.join(output_dir, vcf_filename_noExtension+output_extension)
	vcfHeader = getVCFHeader(vcf_path)
	vcfCols = vcfHeader.split("\t")
	header = '#' + "\t".join(vcfCols[:-1]) + "\t" + "\t".join(tags)
	#write header
	writef(outpath, header)
	# then write data
	tagStrs = []
	for tag in tags:
		tagStrs.append('%INFO/{tag}'.format(tag=tag))
	vcfColStrs = []
	for col in vcfCols:
		vcfColStrs.append('%' + col)
	cmd = "bcftools query -f '{vcfCols}\t{joinedTags}\n' {vcf}".format(vcfCols = "\t".join(vcfColStrs[:-1]), joinedTags="\t".join(tagStrs), vcf=vcf_path)
	print 'vcf to tsv cmd (bcftools extract tags from INFO col): ' + str(cmd)
	subprocess.Popen(cmd, stdout = open_compressed_or_regular(outpath, 'a'), shell=True).wait()
	
	return outpath
	

# generates path for final annotated output file (can be used to find this file if it has already been generated in a previous run)
def generateJoinedOutfilePath(output_dir, sample_name):
	return os.path.join(output_dir, "{sample_name}.annotated.all.tsv".format(sample_name = sample_name))

def tag_exists(input_vcf, tag):
	print 'checking whether tag ' + str(tag) + ' exists in input vcf'
	cmd1 = 'grep "##FORMAT=<ID={tag}" {input_vcf}'.format(tag=tag, input_vcf=input_vcf)
# 	print '1st cmd to check whether tag exists: ' + str(cmd1)
	cmd2 = 'grep "##INFO=<ID={tag}" {input_vcf}'.format(tag=tag, input_vcf=input_vcf)
	proc1 = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	proc2 = subprocess.Popen(cmd2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout1,stderr1 = proc1.communicate()
	rc1 = proc1.returncode
	stdout2, stderr2 = proc2.communicate()
	rc2 = proc2.returncode
	if(rc1 != 0 and rc2 != 0): # not found
		print 'could not find tag ' + str(tag) + ' in input vcf'
		return False
	else:
		print 'found tag ' + str(tag) + ' in input vcf'
		return True
	
	

# extracts desired fields from the FORMAT/INFO column of the input VCF, placing each in a separate column. NOTE: This is different from extracting INFO tags from the dataset VCFs!
def extractFormatCols(inputVCF, out_file_path):
	print 'Extracting FORMAT fields from input VCF'
	# for now, we'll use defaults and worry about YAML spec later
	#Default extract AC, GT, FORMAT/DP, AD (AV and RV are analogous to AD).
	# GT - 1 column per sample (FORMAT/per sample)
	# DP - take mean across samples. 1 column per dataset (FORMAT/per sample)
	# AD - these need to be divided to get the alt allelic fraction. take the mean across samples. 1 column per dataset. 
	# AC - this comes from the INFO tag
	# additional_tag_headers = ['AC', 'GT', 'DP', 'AD']
	# %AC\t%AN
	format_fields = ['AC', 'GT', 'DP', 'AD']
	prefixed_format_fields = []
	for field in format_fields:
		prefixed_format_fields.append('%'+field)
	#debug
	print 'prefixed_format_fields: ' + str(prefixed_format_fields)
	# check that each field exists in vcf before extracting
	actual_prefixed_format_fields = []
	actual_format_fields = []
	for idx,field in enumerate(format_fields):
		if(tag_exists(inputVCF, field)):
			actual_prefixed_format_fields.append(prefixed_format_fields[idx])
			actual_format_fields.append(format_fields[idx])
	print 'actual format fields to query: ' + str(actual_format_fields)
	outfile = open(out_file_path, 'w')
	#print header
	outfile.write("\t".join(actual_format_fields)+"\n")
	outfile.close()
	
	cmd = "bcftools query -f '{formatfields}\n' ".format(formatfields="\t".join(actual_prefixed_format_fields)) + inputVCF # -H prints header but it includes other junk, e.g. # [1]AC	[2]GT	[3]DP	[4]AD
	#+ " | awk '{OFS=\"\t\";print $0,$6/$7}'"
	print 'cmd to extract FORMAT fields from input VCF: ' + cmd
	proc = subprocess.Popen(cmd, shell=True, stdout=open(out_file_path, 'a'), stderr=subprocess.PIPE) # TODO pipe stdout to both python and file (since bcftools gives error as part of stdout, not stderr)
	stderr = proc.communicate()[0]
	if(proc.returncode != 0):
		raise ValueError('Failed to extract format fields from input vcf. Cmd: ' + str(cmd) + '\nError: ' + str(stderr) + '\nReturn code: ' + str(proc.returncode))
	#else
	return out_file_path


# checks that the format cols have appropriate number of lines (matching the input VCF file, including header line and excluding VCF info lines)
def checkFormatCols(inputVCF, format_cols_file):
	return stmp_annotation_checker.check_file_numlines([inputVCF, format_cols_file], ignore_vcf_info_lines=True, raise_exception=True)


def readFilesLines(fileHandlers, delimiter='\t'):
	line_combined = []
	for handler in fileHandlers:
		line = handler.readline().rstrip("\n")
		line_combined.append(line)
	return "\t".join(line_combined)


# joins the files together to form 1 big annotated TSV
# NOTE: requires all files to have exactly the same # of lines in the same order
# Also will not work properly if point or region annotation output files have multiple header lines (which start with #)
def joinFiles(inputVCF, snpeff_tsv, annovar_tsv, region_tsvs, point_tsvs, output_dir, skipJoinCheck, yaml_commands, skip=False, skip_annovar=False, skip_format_fields=True): # by default, don't skip in case the output file was updated	
	print 'Merging all annotations into a single file'
	#debug
	print 'Region annotation dataset to tsv file dict: ' + str(region_tsvs)
	print 'Point annotation dataset to tsv file dict: ' + str(point_tsvs)
	if(not skip_format_fields):
		formatFieldsFile = extractFormatCols(inputVCF, os.path.join(general_utils.get_parent_dir(inputVCF), 'format_fields.tsv'))
		checkFormatCols(inputVCF, formatFieldsFile) # TODO perhaps add user-specifiable option to skip checks
	
	ivh = open_compressed_or_regular(inputVCF, "r")
	if(not skip_format_fields):
		formatFields = open(formatFieldsFile, 'r')
	snpeff = open(snpeff_tsv, "r")
	if(not skip_annovar):
		annovar = open(annovar_tsv, "r")
	
	ordered_dataset_headers = []
	ordered_dataset_paths = []
	datasets_yaml = yaml_utils.get_datasets(yaml_commands)
	for dataset in datasets_yaml:
		dataset_annotated_name = datasets_yaml[dataset][yaml_keys.kDAnnotation]
		if(dataset_annotated_name in point_tsvs):
			dataset_path = point_tsvs[dataset_annotated_name]
		elif(dataset_annotated_name+'_r' in region_tsvs):
			dataset_path = region_tsvs[dataset_annotated_name+'_r']
		else:
			# ignore datasets that we are not importing automatically if they cannot be found
			if(not datasets_yaml[dataset][yaml_keys.kDImportIfMissing]):
				continue
			else:
				raise ValueError('Dataset ' + str(dataset) + ' could not be found in point or region annotated output files. Please ensure that this dataset has been imported into the database.')
		# output this dataset 
		ordered_dataset_headers.append(dataset_annotated_name)
		ordered_dataset_paths.append(dataset_path)
	
	ordered_dataset_handlers = []
	for dataset_path in ordered_dataset_paths:
		ordered_dataset_handlers.append(open_compressed_or_regular(dataset_path, 'r'))
	
	out_filepath = generateJoinedOutfilePath(output_dir, getSampleName(inputVCF))
	output = open(out_filepath, "w")
	
	allCols = {}
	fileCols = []
	foundAnnovarHeader = False
	
	lineNum = -1
	for (vcfLineNum, line) in enumerate(ivh):
		if(line[0:2] != '##'):
			lineNum += 1
			line = line.rstrip("\n")
			
			# remove "." and " "
			
			# add existing cols to hash to prevent duplicate columns in final output
			if(line.startswith('#')):
				for col in line.split("\t"):
					if col.lower() not in allCols:
						allCols[col.lower()] = 1
			
			if(not skip_format_fields):
				formatLine = formatFields.readline()				
			
			snpeffLine = snpeff.readline()
			if(not skip_annovar):
				annoLine = annovar.readline()
			
			if(not skip_format_fields):
				# insert VCF extracted format cols
				line = line + "\t" + formatLine.rstrip("\n")
			
			# check lines to make sure the rows are aligned
			if(not skipJoinCheck):
				# point check
				# TODO update point check code
				print 'Warning: join check currently disabled for point annotations.'
# 				if(not pointLine.startswith('#')):
# 					ptLineStr = pointLine.rstrip("\n").split("\t")[0:7] # chr, coord, id, ref, alt, qual, filter, info
# 					inLineStr = line.split("\t")[0:7]
# 					if(ptLineStr != inLineStr):
# 						print 'Warning joining annotations: input VCF line ' + str(vcfLineNum) + ' does not match point annot line ' + str(lineNum)
# 						print 'point annot line excerpt: ' + str(ptLineStr)
# 						print 'input VCF line excerpt: ' + str(inLineStr)
				
				# annovar check
				if(not skip_annovar):
					if(not annoLine.startswith('#')):
						annoCols = annoLine.rstrip("\n").split("\t")
						chr = annoCols[0]
						start = annoCols[1]
						stop = annoCols[2]
						ref = annoCols[3]
						alt = annoCols[4]
						inCols = line.split("\t")
						if(not (annoCols[0] == inCols[0] # chr
							and inCols[1] == annoCols[1] # start
							#and inCols[1] <= annoCols[2] # stop
							and inCols[3] == annoCols[3] # ref
							and inCols[4] == annoCols[4] # alt
							)
							and annoCols[3] != '-' and annoCols[3] != '.' and annoCols[3] != '0'
							and annoCols[4] != '-' and annoCols[4] != '.' and annoCols[4] != '0'
							and inCols[3] != '-' and inCols[3] != '.' and inCols[3] != '0'
							and inCols[4] != '-' and inCols[4] != '.' and inCols[4] != '0'
						):
							print 'Warning joining annovar functional annotations: input vcf line ' + str(vcfLineNum) + ' does not match annovar line ' + str(lineNum)
							icp = inCols[0:2]
							icp.extend(inCols[3:5])
							acp = annoCols[0:2]
							acp.extend(annoCols[3:5])
							print 'input vcf line (partial): ' + str("\t".join(icp))
							print 'annovar line (partial): ' + str("\t".join(acp))
						
				# region check
				# TODO finish
			
			
			# deal with annoLine separately
			if(not skip_annovar):
				annoLine = annoLine.rstrip("\n")
				if(not foundAnnovarHeader): # first line is header in annovar
					annoLine = annoLine.replace('.', '_')
					annovarHeader = annoLine
					# prepend "Annovar_" to each column in the annovar header
					annovarHeaderCols = annovarHeader.split("\t")
					for idx,col in enumerate(annovarHeaderCols):
						annovarHeaderCols[idx] = 'Annovar_' + col
					annovarHeader = "\t".join(annovarHeaderCols)
					annoLine = annovarHeader # since we print out header from annoLine below
					foundAnnovarHeader = True
				line = line + "\t" + "\t".join(annoLine.split("\t")[5:len(annovarHeader.split("\t"))-1]) # include all columns except the last one (Otherinfo), which spans multiple columns
				
			#deal with snpeff line (just info column for now)
			# TODO extract snpeff annotations instead of outputting entire INFO column
			snpeffLine = snpeffLine.rstrip("\n")
			if(lineNum == 0): # snpeff header
				snpeff_header = snpeffLine.split('\t')
				try:
					snpeff_info_idx = snpeff_header.index('INFO')
				except ValueError as err:
					#attach additional debugging info
					if not err.args:
						err.args=('',)
					err.args += ('snpeff header: ' + str(snpeff_header),)
					err.args += ('snpeff file: ' + str(snpeff_tsv),)
					raise
				line += "\t" + 'snpeff_INFO'
			else: # extract snpeff data
				snpeff_info = snpeffLine.split("\t")[snpeff_info_idx]
				line += "\t" + snpeff_info
				
				
			# deal with point + region annotation files
			for fileIndex,dataset_handler in enumerate(ordered_dataset_handlers):
				dataset_line = dataset_handler.readline().rstrip("\n")
				lineContents = dataset_line.split("\t")
				# deal with other lines
				if(line.startswith('#')): #header - warning: will not work properly if multiple header lines
					outputCols = []
					for idx, col in enumerate(lineContents):
						if(col.lower() not in allCols
						and not col.lower().endswith('chrom') # we already have these
						and not col.lower().endswith('pos')
						and not col.lower().endswith('coord')
						and not col.lower().endswith('ref')
						and not col.lower().endswith('alt')
						and col.lower() != 'chromosome' # same as CHROM
						and col.lower() != 'coordinate' # same as POS
						and col.lower() != 'reference' # same as REF
						and col.lower() != 'alternate' # same as ALT
						and col != ''
						):
							allCols[col.lower()] = 1
							outputCols.append(idx)
							
					fileCols.append(outputCols)	
					
				lineAddArr = []
				
				for i in fileCols[fileIndex]:
					if (i < len(lineContents)):
						lineAddArr.append(lineContents[i])
					else:
						lineAddArr.append('')
				line = line + "\t" + "\t".join(lineAddArr)
			
			#process lines
			line = line.replace('#', '') # '#' messes up import into R
			
			# convert any . to '' in a given entry (to be consistent)
			lineElts = line.split("\t")
			
			# clean up . character in columns (to save space in the final output file)
			for idx, lineElt in enumerate(lineElts):
				if(lineElt.lstrip().rstrip() == '.' or lineElt.lstrip().rstrip() == ''):
					lineElt = ''
					lineElts[idx] = lineElt
			
			## Annotation postprocessing (in-place). Add consensus cols as specified in YAML.
			if(lineNum == 0): #first line = header
				header = list(lineElts) # header DOES NOT include new consensus cols!
				# inject postprocessed consensus col headers
				lineElts = get_consensus_col_header_line(lineElts, yaml_commands)
			else:
				# inject postprocessed consensus cols into annotated output
				lineElts = get_consensus_col_line(lineElts, header, yaml_commands)
			
			line = "\t".join(lineElts)
			output.write(line + "\n")
	
	output.close()
	return out_filepath

# gets the consensus columns in the proper order (look up each one in the YAML for details of what cols it is generated from)
def get_ordered_consensus_cols(consensus_col_cmds):
	consensus_fields = []
	consensus_cmds = consensus_col_cmds
	if(yaml_keys.kAConsensusColumnsOrder in consensus_cmds):
		consensus_fields = consensus_cmds[yaml_keys.kAConsensusColumnsOrder]
	for consensus_field in consensus_cmds: # all other (unordered) consensus fields
		if(consensus_field != yaml_keys.kAConsensusColumnsOrder and consensus_field not in consensus_fields):
			consensus_fields.append(consensus_field)

	return consensus_fields


def get_consensus_col_names(consensus_col_cmds, suffix):
	consensus_col_names = []
	ordered_consensus_cols = get_ordered_consensus_cols(consensus_col_cmds)
	for consensus_field in ordered_consensus_cols:
		consensus_col_name = consensus_field+'_' + suffix
 		consensus_col_names.append(consensus_col_name)
 	return consensus_col_names

# returns header line with consensus cols inserted in appropriate position (according to YAML)
# header = header elements (list), NOT a string!
# same for output (line elements as a list, not a string)
def get_consensus_col_header_line(header, yaml_commands):
	# TODO extract desired VCF FORMAT tags here as well
	# for now, we'll use defaults and worry about YAML spec later
	#Default extract AC, GT, FORMAT/DP, AD (AV and RV are analogous to AD).
	# GT - 1 column per sample (FORMAT/per sample)
	# DP - take mean across samples. 1 column per dataset (FORMAT/per sample)
	# AD - these need to be divided to get the alt allelic fraction. take the mean across samples. 1 column per dataset. 
	# AC - this comes from the INFO tag
	
	annotation_cmds = yaml_commands[yaml_keys.kModules][yaml_keys.kAnnotation]
	consensus_col_cmds = annotation_cmds[yaml_keys.kAConsensusColumns]
	consensus_col_suffix = annotation_cmds[yaml_keys.kAConsensusColumnSuffix]
	offset = annotation_cmds[yaml_keys.kAConsensusColumnOffset]
	
	# insert headers for additional VCF tags
# 	additional_tag_headers = ['AC', 'GT', 'DP', 'AD']
# 	for tag_header in additional_tag_headers:
# 		header.insert(offset, tag_header)
# 		offset += 1
	# now insert the user-specified consensus headers after the VCF tags
	consensus_col_headers = get_consensus_col_names(consensus_col_cmds, consensus_col_suffix)
	for idx,consensus_header in enumerate(consensus_col_headers):
		header.insert(offset+idx, consensus_header)
	
	return header
	
	
# returns line with added consensus cols (in appropriate position as specified in YAML)
# header = header elements (list), NOT a string! index of specific headers must also line up with the lineElts (so give the header BEFORE inserting consensus columns into it)
# same for output (line elements as a list, not a string)
def get_consensus_col_line(lineElts, header, yaml_commands):
	consensus_col_values = []
	annotation_cmds = yaml_commands[yaml_keys.kModules][yaml_keys.kAnnotation]
	consensus_col_cmds = annotation_cmds[yaml_keys.kAConsensusColumns]
	consensus_col_criteria = annotation_cmds[yaml_keys.kAConsensusCriteria]
	offset = annotation_cmds[yaml_keys.kAConsensusColumnOffset]
	# TODO insert VCF tag values before consensus col values
	# TODO finish
# 	vcf_tag_headers = ['AC', 'GT', 'DP', 'AD']
# 	#AC
# 	cmd = "bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC\t%AN\n' " + vcf_file_loc + " | awk '{OFS=\"\t\";print $0,$6/$7}' > " + outfilepath
# 	cmd = 'bcftools query -f %INFO/AC'
	
	# now insert the user-specified consensus column values
	consensus_fields = get_ordered_consensus_cols(consensus_col_cmds)
	# debug
# 	print 'ordered consensus cols: ' + str(consensus_fields)
	for consensus_field in consensus_fields:
		consensus_cols = yaml_utils.convertColumns(consensus_col_cmds[consensus_field], yaml_commands)
		#debug
# 		print 'consensus cols: ' + str(consensus_cols)
		if(consensus_field in consensus_col_criteria): # apply special criteria for finding a consensus beyond default behavior
			consensus_value = get_consensus_value(lineElts, header, consensus_cols, consensus_col_criteria[consensus_field])
		else:
			consensus_value = get_consensus_value(lineElts, header, consensus_cols)
		#debug: print consensus value
# 		print 'consensus value: ' + str(consensus_value)
		consensus_col_values.append(consensus_value)
	
	# now inject the values into the appropriate place in lineElts
	for idx,value in enumerate(consensus_col_values):
		lineElts.insert(offset+idx, str(value))
	
	return lineElts


def get_values(lineElts, header, consensus_cols):
	values = []
	for consensus_col in consensus_cols:
		#debug
# 		print 'header: ' + str(header)
		idx = header.index(consensus_col)
		#debug
# 		print 'idx: ' + str(idx)
		elt = lineElts[idx]
		if ';' in elt:
			elt = elt.split(';')
			values.extend(elt)
# 		elif ',' in elt:
# 			elt = elt.split(',')
# 			values.extend(elt)
		else:
			values.append(elt)
	return values

# for now, we pick the first value from the list that isn't blank or "unknown" (case-insensitive search). If no values match these criteria, we output ''.
def get_consensus_value(lineElts, header, consensus_cols, consensus_criteria='default'):
	values = get_values(lineElts, header, consensus_cols)
	#debug
# 	print 'values: ' + str(values)
	if(consensus_criteria == 'default'):
		for value in values:
			if value != '' and value != ' ' and 'unknown' not in value.lower() and 'NONE' not in value:
				return value
		#else no known value
		return ''
	#else
	if(consensus_criteria == yaml_keys.kAConsensusCriteriaTypeMaxFreq): # get max frequency value (compares all values as floats)
		maxValue = 0.0
		foundHit = False
		for value in values:
			if(value != '' and value != ' '):
				try:
					if(float(value) > float(maxValue)):
						foundHit = True
						maxValue = float(value)
				except ValueError:
					print 'Warning: value ' + str(value) + ' not a float? ' + ' or maxValue ' + str(maxValue) + ' not a float?'
		if(foundHit):
			return maxValue
		else:
			return ''
	#else
	if(isinstance(consensus_criteria, list)): # list of values in order of priority
		valuesDict = general_utils.list2dict(values)
		if(len(valuesDict) == 1 and '' in valuesDict): # all specified columns in that row are empty
			return ''
		for value in consensus_criteria:
			if(value in valuesDict):
				return value
		#else didn't find the value in the list, so fall back to default algo
		print 'warning: could not find any of the specified consensus criteria values ' + str(consensus_criteria) + ' in row values ' + str(values)
		return get_consensus_value(lineElts, header, consensus_cols, 'default') # rerun using default consensus criteria
	#else
	raise ValueError('Unknown consensus criteria specified in YAML: ' + str(consensus_criteria))


def root_or_cwd(mydir):
	# If the user has specified a file using a partial filepath relative to the current working directory, complete that path so that 
	# the file may be located.  If the filepath is absolute (i.e. starting from the root directory) then leave it alone.
	if mydir[0] != '/':
		return os.path.join(os.getcwd(), mydir) # merge the current working directory to the provided partial path
	else:
		return(mydir) # no change


# converts a relative path to an absolute path with respect to the directory where this script is located (NOT the current working directory)
def relativeToAbsolutePath_scriptDir(relativePath):
	script_dir = os.path.dirname(os.path.realpath(__file__))
	if(not relativePath.startswith('/')):
		relativePath = os.path.join(script_dir, relativePath)
	return relativePath

# connect to SQLite database
def connect_db(db_file, host_name='', user_name='', db_name='', unix_socket_loc=''):
	return db_utils.connect_db(db_file, host_name, user_name, db_name, unix_socket_loc)


# gets filename from url
def url2filename(url):
	return url.split('/')[-1]


# converts YAML type (e.g. char) to database type (e.g. varchar(1))
def yamlTypeToDBType(yamlType):
	if(yamlType.lower() == 'character' or yamlType.lower() == 'char'):
		return 'varchar(1)'
	else:
		return yamlType

# looks up SQLite data types for standard VCF columns
def getVCFColTypes(colHeaders):
	colTypes = []
	for header in colHeaders:
		if(header.lower() in vcfHeaders.kVCFColTypes):
			colTypes.append(vcfHeaders.kVCFColTypes[header.lower()])
		else:
			print 'error could not find type of column ' + header
			colTypes.append('')
	
	return colTypes


def is_region_dataset(dataset_name):
	return db_utils.is_region_dataset(dataset_name)

# gets name of column referring to "starting" position (startKey), whether it is called "start", "begin", "pos", etc.
# TODO use this method when forming sql queries instead of inline code to find startKey
def getStartKey(colHeadersLowercase, startKeys):
	for startKey in startKeys:
		if(startKey.lower() in colHeadersLowercase):
			return startKey.lower()
	
	print 'could not get startKey'
	return ''

# whether a given set of column headers has an end position (i.e. refers to a range)
def hasEndKey(colHeaders, endKeys):
	for colHeader in colHeaders:
		if(colHeader.lower() in endKeys):
			return True
	return False

# whether any heaeder in the given list of columns matches a key in the list of keys.
# NOTE: colHeaders are converted to lowercase, keys are NOT (so must be lowercase to begin with).
def hasKey(colHeaders, keys):
	return hasEndKey(colHeaders, keys)

# detects variants of chrom, start, stop, ref, alt and changes them to these standard names for consistency
# NOTE: CURRENTLY UNUSED (standardizeColHeaders below is used instead)
def standardizeYamlHeader(db_yaml):
	colHeadersLc = makeListLowercase(db_yaml['ColumnHeaders'])
	chrIdx = listInLcList(CHR_HEADERS, colHeadersLc)
	if(chrIdx != None):
		db_yaml['ColumnHeaders'][chrIdx] = 'chrom'
	startIdx = listInLcList(START_COLUMN_HEADERS, colHeadersLc)
	if(startIdx != None):
		db_yaml['ColumnHeaders'][startIdx] = 'start'
	stopIdx = listInLcList(STOP_COLUMN_HEADERS, colHeadersLc)
	if(stopIdx != None):
		db_yaml['ColumnHeaders'][stopIdx] = 'stop'
	
	return db_yaml


# detects variants of chrom, start, stop, ref, alt and changes them to these standard names for consistency
def standardizeColHeaders(colHeader):
	colHeadersLc = makeListLowercase(colHeader)
	chrIdx = listInLcList(CHR_HEADERS, colHeadersLc)
	if(chrIdx != None):
		colHeader[chrIdx] = 'chrom'
	startIdx = listInLcList(START_COLUMN_HEADERS, colHeadersLc)
	if(startIdx != None):
		colHeader[startIdx] = 'start'
	stopIdx = listInLcList(STOP_COLUMN_HEADERS, colHeadersLc)
	if(stopIdx != None):
		colHeader[stopIdx] = 'stop'
		
	return colHeader


## special DB import helper functions
# annotates existing clinvar TSV with star ratings
def clinvar_addStars(clinvar_tsv_path, colHeaders, colTypes, yaml_cmds, clinStarHeader=vcfHeaders.kClinvarStarHeader, output_dir=None):
	
	#debug
# 	print 'clinvar tmp tsv path (prior to star annotation): ' + clinvar_tsv_path
	
	filename = getFilename(clinvar_tsv_path)
	filename_noExtension = stripExtension(filename)
	if(output_dir == None):
		output_dir = os.path.dirname(clinvar_tsv_path)
	outpath = os.path.join(output_dir, filename_noExtension+'.stars.tsv')
	infile = open(clinvar_tsv_path, 'r')
	outfile = open(outpath, 'w')
	header = []
	for idx,line in enumerate(infile):
		line = line.rstrip("\n")
		if(idx == 0):
			if(line.startswith('#')):
				line = line[1:]
			header = line.split("\t")
			outfile.write('#' + line + "\t" + clinStarHeader + "\n")
			continue
		#else
		star = vcfUtils.clinvarStars(line, header, yaml_cmds)
		outfile.write(line + "\t" + str(star) + "\n")
	
	colHeaders.append(clinStarHeader)
	colTypes.append('int')
	
# 	outfile.close()
# 	os.chmod(outpath, 'g+rw')
	
	return outpath, colHeaders, colTypes


# loads datasets using info from the YAML file
def setup_db_yaml(database_connection, yaml_commands, skip=True):
	c = database_connection.cursor()
	
	# get existing list of tables
	c.execute('SELECT name FROM sqlite_master WHERE type = "table";')
	tables= c.fetchall()
	table_names = [t[0] for t in tables]
	
	#debug
# 	print 'table names: ' + str(table_names)
	datasets_yaml = yaml_utils.get_datasets(yaml_commands)
	dataset_defaults_yaml = yaml_utils.get_dataset_defaults(yaml_commands)
	for db in datasets_yaml:
		#load default info
		datasetPath = dataset_defaults_yaml[yaml_keys.kDDatasetsPath]
		if(not os.path.isabs(datasetPath)):
			datasetPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), datasetPath)
		if not os.path.exists(datasetPath):
			os.makedirs(datasetPath)
		
		if(db == 'default' or db == yaml_keys.kModules):
			continue
		if(not datasets_yaml[db][yaml_keys.kDImportIfMissing]):
			print 'Skipping ' + str(db) + ' because Import_If_Missing set to False for this dataset in datasets.yml.'
			continue
		
		#else
		startLine = datasets_yaml[db][yaml_keys.kDStartingLine] if yaml_keys.kDStartingLine in datasets_yaml[db] else 0
		# check if dataset already in db and skip if it is
		db_name = datasets_yaml[db][yaml_keys.kDAnnotation].replace('.', '_').replace('-', '_') # get the name of the database from the YAML, replacing illegal characters
		if(datasets_yaml[db][yaml_keys.kDCategory] == yaml_keys.kDCategoryTypeRegion):
			db_name += '_r'
	
		#Skip tables already present (by name, that is) in the database.
		#debug
		print 'db name: ' + db_name
		if(db_name in table_names and skip):
			print 'Skipped completed DB: ' + db_name
			continue
		elif db_name in table_names: # if user has chosen not to skip completed datasets and table name is present in the database
			# drop table
			c.execute("DROP table if exists '{table}';".format(table = db_name))
		
		#handle URL download
		if datasets_yaml[db][yaml_keys.kDSource].lower().startswith('http:') or datasets_yaml[db][yaml_keys.kDSource].lower().startswith('https:') or datasets_yaml[db][yaml_keys.kDSource].lower().startswith('ftp:'):
			filename = url2filename(datasets_yaml[db][yaml_keys.kDSource])
			db_full_path = os.path.join(datasetPath, filename) # Complete the absolute path to the dataset text file.
			if not os.path.isfile(db_full_path):
				db_full_path = downloadDB(datasets_yaml[db], db_full_path, datasetPath, 'URL') # for now, just URL downloads (no Annovar)
		
		#handle file path instead of URL
		else:
			db_full_path = root_or_cwd(datasets_yaml[db][yaml_keys.kDSource])
			if not os.path.isfile(db_full_path):
				db_full_path = os.path.join(datasetPath, datasets_yaml[db][yaml_keys.kDSource])
			if not os.path.isfile(db_full_path):
				raise IOError('Could not find source ' + str(datasets_yaml[db][yaml_keys.kDSource]) + ' in ' + str(root_or_cwd(datasets_yaml[db][yaml_keys.kDSource])) + ' or ' + str(os.path.join(datasetPath, datasets_yaml[db][yaml_keys.kDSource])))
			
			filename = getFilename(db_full_path)
		
		# now do the database import
		c = database_connection.cursor()

		colHeaders = datasets_yaml[db][yaml_keys.kDColumnHeaders] # headers/names of columns in input file
		colTypes = datasets_yaml[db][yaml_keys.kDColumnTypes]
		
		# if file is VCF, convert to TSV with desired INFO tags as columns
		if(filename.endswith('.vcf') or filename.endswith('.vcf.gz')):
			vcfHeader = getVCFHeader(db_full_path)
			# convert to TSV (so each INFO tag goes in a separate column)
			db_full_path = vcf2tsv(db_full_path, colHeaders, datasetPath, output_extension='.tmp.tsv' if datasets_yaml[db][yaml_keys.kDAnnotation]=='clinvar' else '.tsv')
			# update colHeaders and types to reflect standard VCF cols - WARNING: this will not work with multisample vcfs (need to figure out how many sample columns to ignore, not just the last column)
			vcfHeaderCols = vcfHeader.split("\t")[:-1] # ignore INFO column
			vcfTypes = getVCFColTypes(vcfHeaderCols)
			colHeaders = vcfHeaderCols+colHeaders
			colTypes =  vcfTypes+colTypes
			# update filename
			filename = getFilename(db_full_path)
			
		# special extra conversion for clinvar
		if(datasets_yaml[db][yaml_keys.kDAnnotation]=='clinvar'):
			#debug
# 			print 'colheaders: ' + str(colHeaders)
# 			print 'colTypes: ' + str(colTypes)
# 			print 'file: ' + str(db_full_path)
			[db_full_path, colHeaders, colTypes] = clinvar_addStars(db_full_path, colHeaders, colTypes, yaml_cmds=yaml_commands, clinStarHeader=vcfHeaders.kClinvarStarHeader)
			filename = getFilename(db_full_path)
		
		#standardize column headers
		colHeaders = standardizeColHeaders(colHeaders)

		# get dataset column info
		cols = []
		colIndices = {}
		for idx,name in enumerate(colHeaders):
			if(name != ''):
				cols.append("'" + name + "'" + ' ' + "\"{type}\"".format(type=yamlTypeToDBType(colTypes[idx])))
				colIndices[name] = idx
		
		print 'loading ' + db_name
		
		
		colHeaderKeysLower = set(k.lower() for k in colHeaders)
		
		chromKey = 'chrom' if hasKey(colHeaderKeysLower, CHR_HEADERS) else ''
		startKey = 'start' if hasKey(colHeaderKeysLower, START_COLUMN_HEADERS) else ''
		stopKey = 'stop' if hasKey(colHeaderKeysLower, STOP_COLUMN_HEADERS) else ''
		
		# create new table in our db
		c.execute('BEGIN')
		if(is_region_dataset(db_name)):
			c.execute("create table if not exists {db_name}({cols},\
				primary key({chromkey}, {startkey}, {stopkey}));".format(db_name = db_name, stopkey=stopKey, chromkey=chromKey, startkey=startKey, cols = ", ".join(cols)))
			c.execute("create index if not exists start_idx on {db_name}({startkey});".format(db_name = db_name, startkey=startKey))
			c.execute("create index if not exists stop_idx on {db_name}({stopkey});".format(db_name = db_name, stopkey=stopKey))
			c.execute("create index if not exists start_stop_idx on {db_name}({startkey}, {stopkey});".format(db_name = db_name, stopkey=stopKey, startkey=startKey))
		else:
			query = "create table if not exists {db_name}({cols},\
				primary key({chromkey}, {startkey}, ref, alt));".format(db_name = db_name, chromkey=chromKey, startkey=startKey, cols = ", ".join(cols))
			#debug
			print 'db query: ' + query
			c.execute(query)
			c.execute("create index if not exists start_idx on {db_name}({startkey});".format(db_name = db_name, startkey=startKey))

		# load data into our db
		print 'loading file and updating table...'
		db_import(c, db_full_path, db_name, startLine=startLine, colMappings=colIndices, isVCF=filename.endswith('.vcf'), escapeStr='#', enforceTableTypes=False) # header won't be incorrectly imported as long as prefixed with '#'
		database_connection.commit()
		print 'Done importing ' + db_name
		
	# generate BED files for any newly imported region datasets
	c.execute('SELECT name FROM sqlite_master WHERE type = "table";')
	tables= c.fetchall()

	#determine which tables involve ranges (start != stop) and which are single location (start == stop)
 	range_tables = [t[0] for t in tables if t[0].endswith('_r') and not t[0].startswith('sample')]
 	
 	db_bed_dir = relativeToAbsolutePath_scriptDir(yaml_utils.get_dataset_defaults(yaml_commands)[yaml_keys.kDBedPath])
	bed_delimiter = yaml_commands[yaml_keys.kModules][yaml_keys.kAnnotation][yaml_keys.kABedInternalDelimiter]
 	
 	generate_bed_formatted_dbs(database_connection, range_tables, bed_delimiter=bed_delimiter, db_bed_dir=db_bed_dir, yaml_cmds=yaml_commands, force_overwrite=not skip)
 	
 		
# UNUSED and DEPRECATED - replaced by setup_db_yaml
def setup_db(db_dir, database_connection, skip):
	# This script both updates extant datasets and creates new ones within the MySQL database.  
	# A block of code, commented with a line beginning with "format:", accommodates each format
	# in which the incoming datasets may occur.  A dataset is not updated/created when its name
	# appears in the existing tables. 
	c = database_connection.cursor()

	for f in os.listdir(db_dir): # For each text dataset specified in the input directory
		
		# skip all subdirs
		if os.path.isdir(f): 
			print 'skipping subdir ' + f
			continue
		
		if (f.endswith('bcf.gz')):
			# convert to vcf
			db_full_path = os.path.join(db_dir, f) # Complete the absolute path to the dataset text file.
			vcf2allelefreq(db_full_path, db_dir, not skip)
		
		# Note: will the for loop update and get to the newly converted file as well?
		if f.endswith('txt') or f.endswith('tsv') or f.endswith('vcf'): # Using only the .txt and .vcf files...
			db_name = f[0:len(f) - 4].replace('.', '_').replace('-', '_') # Parse the name of the database, replacing illegal characters.
			db_full_path = os.path.join(db_dir, f) # Complete the absolute path to the dataset text file.

			#Skip tables already present (by name, that is) in the database.
			c.execute('SELECT name FROM sqlite_master WHERE type = "table";')
			tables= c.fetchall()
			table_names = [t[0] for t in tables]
			if (db_name in table_names or db_name + '_r' in table_names) and skip: # if the user has chosen to skip completed datasets and the table name is present in the database...
				print("Skipped completed DB: " + db_name)
				continue

			#comprehend format of database.
			with open_compressed_or_regular(db_full_path) as db_txt: # open the database file...
				
				if f.endswith('vcf'):
					print("Updating " + db_name)
					# scan to find header
					line = db_txt.readline()
					print 'Scanning for VCF header' #debug
					while(len(line) < 3 or line[1] == '#' or line [0] != '#'):
						line = db_txt.readline()
						if(len(line) == 0):
							print 'Warning: no VCF header found. Assuming standard columns (chrom, start, id, ref, alt, qual, filter, info)'
							break
					
					#header
					header = line[1:].rstrip().split("\t") # omits leading #
					
					# start with standard columns (the first 8 columns)
					cols = ['chrom varchar(10)', 'start int', 'id varchar(127)', 'ref varchar(512)', 'alt varchar(512)', 'qual varchar(127)', 'filter varchar(127)', 'info varchar(255)']
					colIndices = {'chrom': 0, 'start': 1, 'id': 2, 'ref':3, 'alt':4, 'qual':5, 'filter':6, 'info':7}
					
					# then deal with nonstandard columns
					# note: currently ignores all columns that do not have a corresponding type in the lookup tables
					for index, col in enumerate(header[8:]):
						index += 8 # since we started at col 8 of header
						foundType = False
						if(col.lower() in stmp_consts.col_table_types_absolute): # make the column names lowercase for case-insensitive comparison
							cols.append(col.lower() + ' ' + stmp_consts.col_table_types_absolute[col.lower()])
							colIndices[col.lower()] = index
							foundType = True
						else: # check prefixes dictionary (linear operation)
							for prefix in stmp_consts.col_table_types_prefixes.keys():
								if(col.lower().startswith(prefix)):
									cols.append(col.lower() + ' ' + prefix)
									colIndices[col.lower()] = index
									foundType = True
									break
						
						if not foundType:
							print 'Warning: unknown column ' + col
							
					# load into DB			
					c.execute("BEGIN")
					c.execute("create table if not exists {db_name}({cols},\
					primary key(chrom, start, ref, alt));".format(db_name = db_name, cols = ", ".join(cols)))
					c.execute("create index if not exists start_idx on {db_name}(start);".format(db_name = db_name))
					
					print("loading file and updating table...")
					db_import(c, db_full_path, db_name, 0, colIndices, '#') # start at line 0, using the appropriate column indices, and excluding lines beginning with '#'
					database_connection.commit()
				
				
				elif f.endswith('txt') or f.endswith('tsv'):
					print 'text file: ' + f
					# check for header and process accordingly
					firstLine = db_txt.readline()
					if(firstLine[0] == '#'):
						header = firstLine[1:].rstrip().split("\t")
						
						# annovar popfreq_all
						annovar_popfreq_colHeaders = ['Chr', 'Start', 'End', 'Ref', 'Alt',
													'PopFreqMax', '1000G_ALL', '1000G_AFR',
													'1000G_AMR', '1000G_EAS', '1000G_EUR',
													'1000G_SAS', 'ExAC_ALL', 'ExAC_AFR',
													'ExAC_AMR', 'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE',
													'ExAC_OTH', 'ExAC_SAS', 'ESP6500siv2_ALL', 'ESP6500siv2_AA',
													'ESP6500siv2_EA', 'CG46']
						if(header == annovar_popfreq_colHeaders):
							print("Updating " + db_name)
							c.execute("BEGIN")
							c.execute("create table if not exists {db_name}(chrom varchar(10), start int, stop int, ref varchar(512), alt varchar(512), pop_freq_max float, '1000g_all' float, '1000g_afr' float, '1000g_amr' float, '1000g_eas' float, '1000g_eur' float, '1000g_sas' float, 'exac_all' float, 'exac_afr' float, 'exac_amr' float, 'exac_eas' float, 'exac_fin' float, exac_nfe float, exac_oth float, exac_sas float, esp6500siv2_all float, esp6500siv2_aa float, esp6500siv2_ea float, cg46 float, \
							primary key (chrom, start, stop, ref, alt));".format(db_name = db_name))
							c.execute("create index if not exists start_idx on {db_name}(start);".format(db_name = db_name))
							c.execute("create index if not exists stop_idx on {db_name}(stop);".format(db_name = db_name))
							c.execute("create index if not exists start_stop_idx on {db_name}(start, stop);".format(db_name = db_name))
							print("loading file and updating table...")
							db_import(c, db_full_path, db_name, 1) # start at line 1, ignoring line 0
							database_connection.commit()
						
						else:
							print 'Skipped DB without understood formatting: ' + db_name
							
					else:
						db_txt.seek(0) # go back to the beginning of the db file
						sample = db_txt.readlines(10)[-1].rstrip().split('\t') # ...and sample a few lines from it in order to see its format.
						
						#format: chrom, start, stop, ref, alt, info (where start = stop)
						if len(sample) == 6 and str.isdigit(sample[0]) and str.isdigit(sample[1]) and str.isdigit(sample[2]) and str.isalpha(sample[3]) \
							and str.isalpha(sample[4]) and sample[1] == sample[2]:
							print("Updating " + db_name)
							c.execute("BEGIN")
							c.execute("CREATE TABLE if not exists {db_name}(chrom varchar(10), start int, stop int, ref varchar(512), alt varchar(512), info varchar(255), \
								primary key (chrom, start, stop, ref, alt));".format(db_name = db_name))	# initialize the table if it doesn't already exist
							c.execute("create index if not exists start_idx on {db_name}(start);".format(db_name = db_name))
							c.execute("create index if not exists stop_idx on {db_name}(stop);".format(db_name = db_name))
							c.execute("create index if not exists start_stop_idx on {db_name}(start, stop);".format(db_name = db_name))
							print("loading file and updating table...")
							db_import(c, db_full_path, db_name, 1) # start at line 1, ignoring line 0
							database_connection.commit()
		
						#format: index, chrom with "chr" to be ignored, start, stop, info, ignored1, ignored2, ignored3, ignored4, ignored5
						#e.g. ['590', 'chr1', '760537', '761137', '6_Weak_Enhancer', '0', '.', '760537', '761137', '16776196']
						elif len(sample) == 10 and str.isdigit(sample[0]) and \
							str.isdigit(sample[2]) and str.isdigit(sample[3]) and str.isdigit(sample[4][0]) and str.isdigit(sample[5]) and sample[2] != sample[3]: 
							#label database as containing ranges, rather than single locations.  This will inform the sort of merging which takes place later.
							db_name = db_name + '_r'
		
							print("Updating " + db_name)
							c.execute("BEGIN")
							c.execute("CREATE TABLE if not exists {db_name}(chrom varchar(10), start int, stop int, info varchar(255), \
								primary key (chrom, start, stop, info));".format(db_name = db_name))	# initialize the table if it doesn't already exist
							print("loading file and updating table...")
							db_import(c, db_full_path, db_name, 1, {'chrom':1, 'start':2, 'stop':3, 'info':4})
		
							#remove 'chr'
							print("removing chr character")
							c.execute("UPDATE {table} set chrom = replace(chrom, 'chr', '') where chrom like 'chr%'".format(table = db_name))
							database_connection.commit()
						
						#format: index, chrom with "chr" to be ignored, start, stop, info, info2
						# (e.g. clinvar region file)
						elif len(sample) == 6 and str.isdigit(sample[0]) and str.isdigit(sample[2]) and str.isdigit(sample[3]): # TODO mb check more carefully
							db_name = db_name + '_r'
							print('Updating ' + db_name)
							c.execute("BEGIN")
							c.execute("CREATE TABLE if not exists '{db_name}'(chrom varchar(10), start int, stop int, info varchar(500), \
								primary key (chrom, start, stop, info));".format(db_name = db_name))	# initialize the table if it doesn't already exist
							print("loading file and updating table...")
							db_import(c, db_full_path, db_name, 0, {'chrom':1, 'start':2, 'stop':3, 'info':'4:6'}, enforceTableTypes=False) # ":" = same slicing notation as Python (goes from 4 to 6-1=5 in this case) # WARNING: CANNOT CHECK TYPE (e.g. int) bc chrom contains "chr" (and hence is initially string instead of int)
		
							#remove 'chr'
							print("removing chr character")
							c.execute("UPDATE {table} set chrom = replace(chrom, 'chr', '') where chrom like 'chr%'".format(table = db_name))
							database_connection.commit()
						
						
						# format: chr start stop ref alt sift_score qual? pass?
						elif len(sample) == 8 and str.isdigit(sample[0]) and str.isdigit(sample[1]) and str.isdigit(sample[2]) and str.isalpha(sample[3]) and str.isalpha(sample[4]) and str.isdigit(sample[5]):
# 							db_name = db_name + '_r' # won't work for region annotation bc not VCF (region annotation just uses the info field). Since start and stop appear to always be the same for avsift though, it doesn't matter.
							print 'Updating ' + db_name
							c.execute("BEGIN")
							c.execute("CREATE TABLE if not exists '{db_name}'(chrom varchar(10), start int, stop int, ref varchar(512), alt varchar(512), score double, column7 varchar(10), column8 varchar(10), \
								primary key (chrom, start, stop));".format(db_name = db_name))	# initialize the table if it doesn't already exist
							print("loading file and updating table...")
							db_import(c, db_full_path, db_name, 0, enforceTableTypes=False)
		
							database_connection.commit()
						
				
						# unknown formatting
						else:
							print("Skipped db without understood formatting: " + db_name)
							pass



def makeListLowercase(list1):
	list2 = []
	for item in list1:
		list2.append(item.lower())
	return list2

# checks if ANY of the items in list are in the lc list
# returns index in the lclist or None if not found
def listInLcList(list1, lclist):
	for item in list1:
		if item.lower() in lclist:
			return lclist.index(item.lower())
	return None

# reorders columns so the standard columns (chrom, start, stop, ref, alt) are in the proper order when we create the table in our database
# currently requires STANDARD column names (chrom, start, stop, ref, alt), NOT alternate names (e.g. chr)
def reorderColumns(colNames):
	colNamesLc = makeListLowercase(colNames)
	if('chrom' in colNamesLc):
		idx = colNamesLc.index('chrom')
		colNamesLc.remove('chrom')
		colNamesLc.insert(0, 'chrom')
		elt = colNames.pop(idx)
		colNames.insert(0, elt)
	if('start' in colNamesLc):
		idx = colNamesLc.index('start')
		colNamesLc.remove('start')
		colNamesLc.insert(1, 'start')
		elt = colNames.pop(idx)
		colNames.insert(1, elt)
	if('stop' in colNamesLc):
		idx = colNamesLc.index('stop')
		colNamesLc.remove('stop')
		colNamesLc.insert(2, 'stop')
		elt = colNames.pop(idx)
		colNames.insert(2, elt)
	if('ref' in colNamesLc):
		idx = colNamesLc.index('ref')
		colNamesLc.remove('ref')
		colNamesLc.insert(3, 'ref')
		elt = colNames.pop(idx)
		colNames.insert(3, elt)
	if('alt' in colNamesLc):
		idx = colNamesLc.index('alt')
		colNamesLc.remove('alt')
		colNamesLc.insert(4, 'alt')
		elt = colNames.pop(idx)
		colNames.insert(4, elt)
		
	return colNames


# generates bed files from region tables in our database so we can use bedtools to perform the region annotation
def generate_bed_formatted_dbs(database_connection, table_names, db_bed_dir, bed_delimiter, yaml_cmds, force_overwrite=False):
	print("Generating bed-formatted tables for region annotation")
	
	db_bed_dir = relativeToAbsolutePath_scriptDir(db_bed_dir)
	
	if(not os.path.exists(db_bed_dir)):
		os.makedirs(db_bed_dir)	

	c = database_connection.cursor()
	
	tables_and_files = [(table, os.path.join(db_bed_dir, table + '.bed')) for table in table_names]
	
	processes = []
	
	for t,f in tables_and_files:
		if not os.path.isfile(f) or force_overwrite:
			print("\nGenerating BED for " + t)
			columnsInfo = getColumnsOfTable(c, t)
			colNames = []
			# get column nmes
			for columnInfo in columnsInfo:
				colNames.append(str(columnInfo[1]))
			# check for standard columns (chrom, start, stop, ref, alt) and add them in proper order
			colNames = reorderColumns(colNames)
			
			# for bedtools, only pull the rows that have genomic coordinates and positive start coordinate
			cmd = "SELECT {colnames} from '{t}' where start != '' and stop != '' and start >= 0;".format(t=t, f=f, colnames = ', '.join(colNames))
			print 'BED generation sql cmd: ' + str(cmd)
			
			# write header - NOTE: DOES NOT INCLUDE CHR, START, STOP (only the columns after these, aka the cols that need to be split up when parsing the intersected file)
			headerf = getHeaderFile(f)
			headerout = open(headerf, 'w')
			headerout.write("\t".join([t+'_'+colname for colname in colNames[3:]]) + "\n")
			
			unsortedf = f+'.unmerged.bed'
			with open(unsortedf, 'w') as outfile:
				for ridx,row in enumerate(c.execute(cmd)):
					strings = []
					infoArr = []
					#debug
					if(ridx < 5):
						print 'row: ' + str(row)
					for idx,x in enumerate(row):
						if(idx == len(row) -1 and (type(x) is str or type(x) is unicode)):
							#debug
# 							print 'x before rstrip: ' + str(x)
							x = x.rstrip("\n\n") # effectively the same as rstrip("\n")
							#debug
# 							print 'x after rstrip: ' + str(x)
						if type(x) is int or type(x) is float or type(x) is str:
							x2 = str(x)
							x2 = x2.encode('utf8')
						else:
							x2 = x.encode('utf8')
						if(idx < 3): # chr, start, stop
							strings.append(x2)
						else: # extra info -> we'll put it all in the "name" column of the BED file
							infoArr.append(x2)
							
					#make name column
					name = bed_delimiter.join(infoArr)
					strings.append(name)
					
					#now write our 4 columns to BED file
					outstr = "\t".join(strings)
					outfile.write(outstr + "\n")
				outfile.close()
			# merge regions during loading so we don't need to use condense_intersectBed_output.py later
			#debug
			print 'bedfile: ' + str(unsortedf)
			
			# SORT AND MERGE WORKAROUND (should work on Mac OS, but does not actually sort or merge the output BED files). WARNING: This will NOT behave as expected unless the condense_intersectBed script is used as part of range annotation
# 			cmd = 'cat "{bedfile}"'.format(bedfile=unsortedf)
# 			print 'bed ignore sort command: ' + str(cmd)
# 			subprocess.Popen(cmd, stdout=open(f, 'w'), shell=True).wait()
# 			print 'unsorted output bedfile: ' + str(f)
			
			# REGULAR SORT AND MERGE (doesn't work on Mac OS bc GNU sort accepts different parameters)
			cmd = 'sort -k1,1V -k2,2n "{bedfile}"'.format(bedfile=unsortedf)
			print 'bed sort command: ' + str(cmd)
			subprocess.Popen(cmd, stdout=open(f, 'w'), shell=True).wait()
			print 'sorted output bedfile: ' + str(f)
# 			print 'sorted bedfile: ' + str(f+'.sorted.bed')
# 			cmd = 'bedtools merge -c 4 -o collapse -delim "{bed_multimatch_internal_delimiter}"'.format(bed_multimatch_internal_delimiter=yaml_cmds[yaml_keys.kModules][yaml_keys.kAnnotation][yaml_keys.kABedMultimatchInternalDelimiter])
# 			print 'bed region merge cmd: ' + str(cmd)
# 			subprocess.Popen(cmd, stdin=open(f+'.sorted.bed', 'r'), stdout=open(f, 'w'), shell=True).wait()
# 			print 'output bedfile: ' + str(f)
# 			processes.append(subprocess.Popen(cmd, stdin=open(f+'.sorted.bed', 'r'), stdout=open(f, 'w'), shell=True))
			# TODO error check bedtools merge and other output for each third-party command (return codes and stderr)
			
		else:
			print("Reusing existing BED for " + t)
		
		# TODO check whether final output BED file has size > 0
			
	database_connection.commit()
	for process in processes:
		print 'bed region merging: waiting on ' + str(process)		
	print 'Done generating BED files.'

# exception handling
def PrintException():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print 'EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj)


# gets our generated header file for the given bed file (file with the same name as the BED file but ends in .bed.header instead of .bed)
def getHeaderFile(bed_file_path):
	return bed_file_path.replace('.bed', '.bed.header')


def expandBED(bedPath, bed_multimatch_internal_delimiter, bed_delimiter, dataset_multimatch_delimiter, fileHeaders):
	print 'Expanding file ' + str(bedPath)
	outSuffix='.tsv'
	f = open(bedPath, 'r')
	out_file_path = bedPath+outSuffix
	out_file = open(out_file_path, 'w')
	prevTime = datetime.datetime.now().second
	for currLine,fl in enumerate(f):
		prevTime = print_progress_timed('line ' + str(currLine), prevTime)
		lineContents = []
		fl = fl.rstrip("\n")
		if(fl == '.' or len(fl) == 0):
			# extend lineContents by appropriate # of blank columns
			lineContents.extend(['']*len(fileHeaders[bedPath]))
		elif(bed_multimatch_internal_delimiter in fl):
			# split annotations and then recombine using specified default delimiter
			annotations = fl.split(bed_multimatch_internal_delimiter)
			finalCols = ['']*len(fileHeaders[bedPath]) # ensures the appropriate number of cols for future hits (which will be appended to the appropriate empty string)
			for annotation in annotations:
				annotationComponents = annotation.split(bed_delimiter)
				if(len(annotationComponents) > 0):
					for idx,component in enumerate(annotationComponents):
						try:
							if(finalCols[idx] != ''):
								finalCols[idx] += dataset_multimatch_delimiter
							finalCols[idx] += component
						except IndexError:
							print 'error: index ' + str(idx) + ' beyond size of finalCols (' + str(len(finalCols)) + ') on file ' + str(bedPath) + ' line ' + str(currLine)
							raise 
							
			lineContents.extend(finalCols)
		else:
			flc = fl.split(bed_delimiter)
			lineContents.extend(flc)
		out_file.write("\t".join(lineContents) + "\n")
	return out_file_path
	
# returns dictionary of annotated dataset name -> path to annotated TSV file.
def annotate_range(database_connection, vcf_file_loc, output_dir_loc, modules_yaml_path, datasets_yaml_path, yaml_commands, skip=False, print_range_cmd_only=False, presorted_vcf=False, force_overwrite_beds=False):
	vcf_numSampleCols = vcfUtils.getNumSampleCols(vcf_file_loc)
	intersected_tmp_dir = os.path.join(output_dir_loc, 'intersected_beds') # TODO make scratch dir a separate param for this function
	db_bed_dir = relativeToAbsolutePath_scriptDir(yaml_utils.get_dataset_defaults(yaml_commands)[yaml_keys.kDBedPath])
	bed_delimiter = yaml_commands[yaml_keys.kModules][yaml_keys.kAnnotation][yaml_keys.kABedInternalDelimiter]
	bed_multimatch_internal_delimiter = yaml_commands[yaml_keys.kModules][yaml_keys.kAnnotation][yaml_keys.kABedMultimatchInternalDelimiter]
	
	c = database_connection.cursor()
	sample_name = 'sample_' + vcf_file_loc[0:len(vcf_file_loc) - 4].replace('.', '_').replace('-', '_').split('/')[-1]
	out_file_loc = os.path.join(output_dir_loc, sample_name + "_annotated.tsv")
	
	print("Annotating variants for sample " + sample_name)

	c.execute('SELECT name FROM sqlite_master WHERE type = "table";')
	tables= c.fetchall()

	#determine which tables involve ranges (start != stop) and which are single location (start == stop)
	# TODO use datasets.yml for range tables instead of database (need to link back to datasets.yml names)
 	range_tables = [t[0] for t in tables if t[0].endswith('_r') and not t[0].startswith('sample')]

 	#perform BEDtools intersections with the range tables
 	print("Performing BEDtools intersections with {tabs} tables".format(tabs = len(range_tables)))

 	processes = []
 	annotated_dataset_names = []
#  	yaml_dataset_names = []
 	out_file_paths = []
	
	if not os.path.exists(intersected_tmp_dir):
		os.makedirs(intersected_tmp_dir)

 	for t in range_tables:	
		out_file_loc = os.path.join(intersected_tmp_dir, sample_name + '_isect_' + t + '.tab')
		bed_file_path = os.path.join(db_bed_dir, t + '.bed')
#  		yaml_dataset_name = yaml_utils.annotated_to_yaml_dataset_name(t, yaml_commands)
 		
 		try:
 			# TODO have condense_intersectbed do the cutting to cut out the VCF cols. That would allow us to easily include multiple cols from the intersectbed output (as opposed to just cutting 1).
			cmd = "(cat {header}; intersectBed -a {vcf_file_loc} -b {fname} -loj | {condense} | cut -f{colIndex})".format(
				colIndex=13+vcf_numSampleCols,
				header = getHeaderFile(bed_file_path),
				vcf_file_loc = vcf_file_loc, 
				fname = bed_file_path,
				condense = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'condense_intersectBed_output.py --modules={modules_yaml}'.format(modules_yaml=modules_yaml_path))
			)
	 		
	 		print 'range annotation cmd: ' + cmd
			processes.append(subprocess.Popen(cmd, stdout = open(out_file_loc, 'w'), shell = True))
# 			yaml_dataset_names.append(yaml_dataset_name)
			annotated_dataset_names.append(t)
			out_file_paths.append(out_file_loc)
			
		except:
			PrintException()
	
	if(print_range_cmd_only):
		return
	
	prevTime = datetime.datetime.now().second
	for p in processes:
		prevTime = print_progress_timed('Waiting on intersectbed range annotation process ' + str(p), prevTime)
		p.wait()
	
	# now merge the range annotations
	annotations = out_file_paths
	
	# check line counts
	stmp_annotation_checker.check_file_numlines(annotations) # TODO add VCF linecount (non ## lines but include header line) to compare against
	
	out_filepath = os.path.join(output_dir_loc, "{name}.range_joined.vcf".format(name = sample_name))
	
	# for now, just using the default delimiter. Later, check if there's a specific delimiter for a given dataset and use that instead.
	dataset_multimatch_delimiter = yaml_utils.get_dataset_defaults(yaml_commands)[yaml_keys.kDMultimatchDelimiter]
	
	
	# expands bed files to have proper # of columns for easy merge into the final annotated file
	# NOTE: expanding .tab files, NOT .bed files!
	def expandBEDs(bedPaths):
		script_dir = general_utils.get_code_dir_abs()
		expand_script_path = os.path.join(script_dir, 'expand_intersectedBeds.py')
		processes = []
		bedOutFiles = []
		outSuffix = '.tsv'
		for bedPath in bedPaths:
			cmd = 'python "{expandScript}" --tab_file "{bedPath}" --out_suffix "{outSuffix}" --config_modules "{config_modules}" --config_datasets "{config_datasets}"'.format(expandScript = expand_script_path, bedPath=bedPath, outSuffix=outSuffix, config_modules=modules_yaml_path, config_datasets=datasets_yaml_path)
			print 'expand intersectedBeds cmd: ' + str(cmd)
			processes.append(subprocess.Popen(cmd, shell=True))
		for idx,process in enumerate(processes):
			print 'waiting on expand beds cmd for ' + str(bedPaths[idx])
			stderr = process.communicate()[0] # check
			returncode = process.returncode
			if(returncode != 0):
				raise ValueError('expand beds cmd failed (return code ' + str(returncode) + ')') #+ str(stderr))
			bedOutFiles.append(bedPaths[idx]+outSuffix)
		return bedOutFiles
	
	# UNUSED
	def mergeFiles(expandedBedPaths, outFilePath):
		cmd = 'paste {expandedBedFiles}'.format(expandedBedFiles=' '.join(expandedBedPaths))
		print 'range annotation mergeFiles cmd: ' + str(cmd)
		print 'writing cmd output to ' + str(outFilePath)
		process = subprocess.Popen(cmd, shell=True, stdout=open(outFilePath, 'w'), stderr = subprocess.PIPE)
		err = process.communicate()[0] # check
		returncode = process.returncode
		print 'cmd return code: ' + str(returncode) + ' stderr: ' + str(err)
		if(returncode != 0):
			raise ValueError('Merging range annotation files (paste command) failed.\nReturn code: ' + str(returncode) + '\nErr: ' + str(err) + '\nCmd: ' + str(cmd))
		# TODO might need to close output file handle here?
	
	print 'merging range annotation .tab files'
	expandedBeds = expandBEDs(annotations)
	# check line counts
	stmp_annotation_checker.check_file_numlines(expandedBeds)
# 	mergeFiles(expandedBeds, out_filepath)
	
	print 'Finished bedtools region (range) annotation'
	
	# generate dictionary of dataset name -> expanded intersected bed path for each intersected file.
	dataset_path_dict = {}
	for idx,annotated_dataset_name in enumerate(annotated_dataset_names):
		dataset_path_dict[annotated_dataset_name] = expandedBeds[idx]
	
	return dataset_path_dict
# 	return expandedBeds
# 	return out_filepath


# escapes the elts in the array with single quotes ''
def escapeElts(array):
	array2 = []
	for elt in array:
		array2.append("'" + elt + "'")
	
	return array2


# prints the given string every [interval] s or longer (depending on how often this method is called). Call it each time the progress (e.g. line number) changes.
# NOTE: prevTime and interval must be in SECONDS
def print_progress_timed(progressStr, prevTime, interval=5):
	currtime = datetime.datetime.now().second
	if(abs(currtime - prevTime) > interval):
		print str(datetime.datetime.now()) + ': ' + progressStr
		return currtime
	#else
	return prevTime

#returns an array with info about each column (column[1] = name, column[2] = type)
def getColumnsOfTable(db_connection_cursor, tablename):
	return db_utils.getColumnsOfTable(db_connection_cursor, tablename)

# Helper method to import the contents of a TSV/VCF/vcf.gz file to a given table in the database
# colMappings = hash mapping field names in the table to column #s in the file (starting at 0)
# e.g. "chrom" -> 1
# escapeStr = escape string for input data file (so it can ignore these lines) -> e.g. '#'
# printDebug = whether or not to print debug messages
# TODO: add option to do something other than replace if the elt to be inserted already exists
def db_import(db_connection_cursor, filename, tablename, startLine=0, colMappings={}, escapeStr='', isVCF=False, printDebug=False, enforceTableTypes=True, stripChrIfPresent=True):
# 	if(printDebug):
	#debug
	print 'importing filename: ' + filename
		
	db_connection_cursor.execute('pragma table_info({table});'.format(table=tablename))
	columns = db_connection_cursor.fetchall();
	colTypes = []
	for column in columns:
		name = column[1]
		colType = column[2]
		colTypes.append(colType)
		
	entries = [] # will be list of tuples to add to db
	
	colNames = colMappings.keys()
	
	dumpedEntries = False

	# imports all current entries into the db, clearing the entries array afterward
	# uses colNames from outer function scope
	def importEntries(entries):	
		colNameStr = ''
		if(len(colNames) > 0):
			colNameStr = '(' + ", ".join(escapeElts(colNames)) + ")"
		db_connection_cursor.executemany("insert or replace into {tablename}{colNameStr} values ({vals});"
										.format(tablename=tablename, 
											colNameStr=colNameStr,
											vals=", ".join("?" for i in range(0, len(entries[0])))), 
										entries)
		entries = [] # clear entries array since done importing
		return entries

	with open_compressed_or_regular(filename, 'r') as f:
		for i in range (0, startLine):
			next(f) # skip lines until we reach the desired start point
		
		idx = 0
		prevTime = datetime.datetime.now().second
		
		for line in f:
			#strip trailing newline
			line = line.rstrip("\n")
			
			# print our progress (line #)
			prevTime = print_progress_timed('importing line ' + str(idx) + ' of ' + tablename, prevTime)
			idx += 1;
			
			line = line.decode('utf-8', 'ignore') # in case there are any byte-strings, convert them to utf-8 or ignore them if they can't be converted. This avoids potential errors when trying to add bytestrings to the sqlite db.
			if(len(escapeStr) == 0 or line[0:len(escapeStr)] != escapeStr):
				lineContents = line.split("\t")
				if printDebug:
 					print 'got line: ' + line 					
				if(len(colNames) > 0): # extract specific columns/VCF tags from file
					colVals = []
					for colName in colNames:
						if(colName.lower() in CHR_HEADERS):
							#strip 'chr' prefix if present in line data (e.g. 'chr1' becomes just '1')
							if(lineContents[colMappings[colName]].startswith('chr')):
								lineContents[colMappings[colName]] = lineContents[colMappings[colName]][3:]
						if(isinstance(colMappings[colName], (int))):
							colVals.append(lineContents[colMappings[colName]]) # add column values in appropriate order
							
						else: # colMappings[colName] is of the form #:# to indicate a range of columns that should be merged using | to separate them
							col_range = colMappings[colName]
							col_range = col_range.split(":")	
							colVals.append("|".join(lineContents[int(col_range[0]):int(col_range[1])])) # TODO maybe make the delimiter customizable			
					
					# verify that columns are of valid type, if type enforcing enabled
					addtoEntries = True
					
					if(enforceTableTypes):
						for idx,colVal in enumerate(colVals):
							if(
							(colTypes[idx].startswith('int') and not colVal.isdigit())
							):
								addtoEntries = False
								break
					
					if(addtoEntries):
						colTuple = tuple(colVals)
						entries.append(colTuple)
						if printDebug:
							print 'appended tuple ' + str(colTuple) + 'to entries'
							
						# bunch into chunks of 1000 entries for improved load performance
						if(len(entries) > 1000):
							#debug
							if not dumpedEntries:
# 								print 'entries: ' + str(entries)
# 								print 'importing ' + str(len(entries)) + ' entries'
								dumpedEntries = True
							entries = importEntries(entries)
				
				else: # add all columns from file
					addtoEntries = True
					
					if(enforceTableTypes):
						for idx,colVal in enumerate(lineContents):
							if(
							(colTypes[idx].startswith('int') and not colVal.isdigit())
							):
								print 'warning: ' + str(colVal) + ' is not of type ' + str(colTypes[idx]) + ' and will be ignored for row ' + str(lineContents)
								addtoEntries = False
								break
							
					if(addtoEntries):
						colTuple = tuple(lineContents)
						entries.append(colTuple)
						if(len(entries) > 1000):
							#debug
# 							print 'importing entries ' + str(entries)
							entries = importEntries(entries)
							
			
			else: # line has escape str
				None
	
	# import any remaining entries (<=1000)
	if(len(entries) > 0):
		#debug
# 		print 'importing entries ' + str(entries)
		entries = importEntries(entries)


#drops a single sample table (used for dropping the sample table after annotation, if specified)
def drop_sample(database_connection, vcf_file_loc):
	sampleTableName = getSampleTableName(vcf_file_loc)
	print 'Dropping sample table ' + sampleTableName + ' from the database'
	c = database_connection.cursor()
	c.execute("DROP table if exists '{table}';".format(table = sampleTableName))

# drops all tables beginning with "sample_"
def drop_samples(database_connection):
	print 'Dropping all sample tables from database'
	
	c = database_connection.cursor()
	c.execute('SELECT name FROM sqlite_master WHERE type = "table";')
	tables= c.fetchall()
	for table in tables:
		if(table[0].startswith('sample')):
			c.execute('drop table {table};'.format(table=table[0]))
	

# perform point annotation
# returns dictionary of annotated dataset name -> path to annotated output TSV file.
def annotate_point(database_connection, vcf_file_loc, output_dir_loc, sample_db_path, print_sql_query_only=False, debug=False):
	c = database_connection.cursor()
	sample_name = 'sample_' + vcf_file_loc[0:len(vcf_file_loc) - 4].replace('.', '_').replace('-', '_').split('/')[-1]
	out_file_loc = os.path.join(output_dir_loc, sample_name + ".point_annotated.tsv")

	print("Annotating variants for sample " + sample_name)
	
	# attach sample db
	c.execute('attach "{sample_db}" as sample;'.format(sample_db=sample_db_path))

	c.execute('SELECT name FROM sqlite_master WHERE type = "table";')
	tables= c.fetchall()

	# perform SQL join with locus tables
	locus_tables = [t[0] for t in tables if not t[0].endswith('_r') and not t[0].startswith('sample')]
	print("Performing SQL join with {tabs} tables:".format(tabs = len(locus_tables)))
	print('\n'.join(locus_tables))

	cmd = '' # initialize command text
	
	#returns whether or not to include the given table (reused in 3 different places)
	def include_table(table):
		return not table.startswith('sample')
	
	sql_queries = []
	tables_for_queries = []
	output_files_by_table = {}
	
	for table in locus_tables:
		if include_table(table):
			tables_for_queries.append(table)
			cmd = ''
			columns = db_utils.getColumnNamesOfTable(c, table)
			for column in columns:
				column_lower = column.lower()
				if(column_lower != 'chrom' and column_lower != 'start' and column_lower != 'stop'
				and column_lower != 'ref' and column_lower != 'alt'
				and column_lower != 'id' and column_lower != 'qual'
				and column_lower != 'filter'
				):
					if(cmd == ''):
						cmd = "select '{db}_{column}'".format(db = table, column=column)
					else:
						cmd += ", '{db}_{column}'".format(db = table, column=column)
			cmd += ' UNION ALL'
			# select output columns
			first = True
			for column in columns:
				column_lower = column.lower()
				if(column_lower != 'chrom' and column_lower != 'start' and column_lower != 'stop'
				and column_lower != 'ref' and column_lower != 'alt'
				and column_lower != 'id' and column_lower != 'qual'
				and column_lower != 'filter'
				):
					#regular
					if(first):
						cmd += " select ifnull({table}.'{column}','') as '{table}_{column}'".format(table = table, column=column) # add tab to keep columns separated
						first = False
					else:
						cmd += ", ifnull({table}.'{column}','') as '{table}_{column}'".format(table = table, column=column) # add tab to keep columns separated
			#perform join
			cmd += ' from ' + sample_name
			cmd += " left join {next_table} on {sample}.chrom = {next_table}.chrom and {sample}.start = {next_table}.start and {sample}.ref = {next_table}.ref and {sample}.alt = {next_table}.alt;".format(next_table = table, sample = sample_name)
			
			sql_queries.append(cmd)
	
	if(print_sql_query_only):
		print 'Point annotation SQL queries:\n' + '\n'.join(sql_queries)
		return
	#else, execute these queries
	for idx,cmd in enumerate(sql_queries):
		print 'executing point annotation sql query: ' + cmd
		table = tables_for_queries[idx]
		out_file_loc = os.path.join(output_dir_loc, table + '.point_annotated.tsv')
		output_files_by_table[table] = out_file_loc
		# remove output file if it already exists
		if os.path.isfile(out_file_loc):
			call(["rm", out_file_loc])
		print 'writing point annotation output file ' + out_file_loc
		with open(out_file_loc, 'w') as f:
			if(idx == 0):
				timer = datetime.datetime.now().second
				rowNum = 0
			for row in c.execute(cmd): # seems to take a while on some input files
				if(idx == 0): # print progress for first output file only
					currtime = datetime.datetime.now().second
					if(abs(currtime - timer) > 5):
						print str(datetime.datetime.now()) + ': writing line ' + str(rowNum) + ' to first point annotation file'
						timer = currtime
				# write each row of the query result to the file (tab-delimited) 
				try:
					f.write(("\t".join(x.encode('utf8').strip("\n") if isinstance(x, basestring) else str(x) for x in row) + "\n"))
				except UnicodeEncodeError:
					print 'UnicodeEncode error in row: ' + str(row)
					raise
				except UnicodeDecodeError:
					print 'UnicodeDecode error in row: ' + str(row)
					raise
				rowNum += 1
	
	print 'Finished writing point annotation output files at ' + str(output_dir_loc)
	
	# return a list (dictionary) of output files by table name
	return output_files_by_table


def get_sample_db_path(db_dir, sample_name):
	return os.path.join(db_dir, sample_name+'.sqlite')

def upload_vcf(database_connection2, vcf_file_loc, db_dir, force_overwrite=False):
	sample_name = getSampleTableName(vcf_file_loc)
	sample_db_path = get_sample_db_path(db_dir, sample_name)
	database_connection = sqlite3.connect(sample_db_path)
	c = database_connection.cursor()

	c.execute('SELECT name FROM sqlite_master WHERE type = "table";')
	tables= [t[0] for t in c.fetchall()]

	if sample_name not in tables or force_overwrite:
		print("Uploading VCF to database for sample " + sample_name)
		c.execute("BEGIN")
		c.execute("DROP table if exists '{table}';".format(table = sample_name))
		c.execute("CREATE TABLE '{table}'( \
			chrom varchar(10), \
			start int, \
			id varchar(127), \
			ref varchar(512), \
			alt varchar(512), \
			qual varchar(127), \
			filter varchar(127), \
			info varchar(255), \
			primary key(chrom, start, ref, alt));".format(table = sample_name))	# initialize the table if it doesn't already exist

		print("updating VCF table...")
		db_import(c, vcf_file_loc, sample_name, 1, {'chrom':0, 'start':1, 'id':2, 'ref':3, 'alt':4, 'qual':5, 'filter':6, 'info':7}, '#')
		database_connection.commit()
	else:
		print(sample_name + " found in database.  Beginning annotation without re-upload.")
	
	return sample_db_path


