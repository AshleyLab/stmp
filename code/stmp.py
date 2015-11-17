#!/usr/bin/env python

'''
Author: Eli Moss elimoss@stanford.edu, Prag Batra prag@stanford.edu

Purpose: 

	a) Provide facile creation/maintenance of SQLite database containing variant annotation data.
	b) Annotate variants (as VCF) with annotations contained in database.
	c) Prioritize annotated variants into 4 different tiers according to functional and other annotations.

Explanation:

	Organization and maintenance of annotation data is complex when many dozens of datasets are involved.  Subsequently, applying
	those annotations to a collection of variants becomes an intricate task rapidly.  This script solves both problems by creating
	and utilizing a SQLite annotation database.  Further detail is provided in the code.

Example:

stmp.py --update_db --config=config/datasets.yml # loads all datasets specified in the YAML file into our database (if the config parameter is not specified it defaults to the above)

stmp.py --vcf={input vcf or vcf.gz file} --output_dir={output directory} # applies annotation and tiering to the input VCF, outputting it in the specified output directory. Uses the yaml files in config/ for annotation. If you wish to use different configuration files, you can either modify these files or copy them to a different location and link to them via the --config and --modules parameters (config should point to datasets.yml and modules to modules.yml).

For a complete list of parameters, just run stmp.py with no options.


'''


import argparse
import stmp_annotation_util
import stmp_tiering_util
import yaml_utils
import pgxUtils
import general_utils
import stmp_annotation_checker # for checking annotated output file

import multiprocessing
from multiprocessing import pool
import time
from time import strftime, localtime
import os
import sys
import yaml
import yaml_keys
import subprocess
import datetime
import gzip

cur_path = os.path.dirname(os.path.abspath(__file__))
resources_path= os.path.join(cur_path, os.path.pardir, 'data')

# parse commandline options.
parser = argparse.ArgumentParser(description = 'Sequence to Medical Phenotypes')

parser.add_argument("--config", dest='yaml', help = 'yaml-formatted annotation/db update specification', default=os.path.join(cur_path, 'config', 'datasets.yml'))
parser.add_argument('--modules', help='yaml-formatted options for each module (annotation, tiering, pgx, trio)', default=os.path.join(cur_path, 'config', 'modules.yml'))

# arguments related to updating the database
db_update_group = parser.add_argument_group("Database initialization/update", 'Initialize DB, or update if already extant.')
db_update_group.add_argument("--update_db", action='store_true', 
	help = 'Initialize or update MySQL database from text annotation files.  Please provide either YAML configuration file or directory containing annotation database text files (see below parameters)')
db_update_group.add_argument('--input_directory', help='Directory containing annotation database text files.')
db_update_group.add_argument('--force_overwrite', action='store_true', help='Forces overwrite of all existing datasets.')
db_update_group.add_argument("--download_datasets_only", action='store_true', help='just download the datasets specified in the YAML file without loading them into the DB or performing annotation')
db_update_group.add_argument('--check_datasets_only', action='store_true', help='just check which datasets from the YAML were downloaded already')
db_update_group.add_argument('--dataset_output_dir', help='output dir for the downloaded datasets')

# arguments related to annotating a vcf
annotate_group = parser.add_argument_group('Annotation Parameters', 'Arguments required in order to perform VCF annotation.')
annotate_group.add_argument("--vcf", help = 'input VCF or gzipped VCF file')
annotate_group.add_argument("--output_dir", help = 'location of output files')
annotate_group.add_argument("--clean_beds" , action='store_true', dest='clean_beds', help="do a clean run instead of reusing previously generated bed files for range annotation (note that generating these files takes a while)")
annotate_group.add_argument("--skip_multiallelic_fix", action='store_true', dest='skip_multiallelic', help="don't use bcftools to split lines with multiple allele changes")
annotate_group.add_argument("--reuse_multiallelic", action='store_true', help='reuse multiallelic split file if it already exists')
annotate_group.add_argument("--force_input_overwrite", action='store_true', dest='force_input_overwrite', help='forces input VCF to be reimported into the DB, even if it has already been imported in the past')
annotate_group.add_argument('--log', help='name of logfile')
annotate_group.add_argument('--database_file', help='path to SQLite database file (default: ../db/annotationDB.sqlite)', default=os.path.join(cur_path, '..', 'db' ,'annotationDB.sqlite'))
annotate_group.add_argument('--annotate_only', action='store_true', help='just do annotation (no tiering, pgx, etc.)')
annotate_group.add_argument('--drop_samples', action='store_true', help='drops sample tables (use this before distributing the database)')
annotate_group.add_argument('--drop_sample', action='store_true', help='drops this sample from the database after completing annotation')
annotate_group.add_argument('--region_annotations_only', action='store_true', help='just do region annotation')
annotate_group.add_argument('--print_region_cmd_only', action='store_true', help="just print region annotation command (don't actually do the annotation or anything else")

# variant tiering (prioritization)
tiering_group = parser.add_argument_group('Variant Prioritization', 'Arguments for variant prioritization')
tiering_group.add_argument('-t', '--target_genes', default=None, help='target candidate gene list (1 gene per line, if TSV make sure the gene name is in the first column) for filtering variants prior to prioritizing them.')
tiering_group.add_argument('-e', '--ethnicity', default =None, help = 'ethnicity for allele frequency filtering')
tiering_group.add_argument('-s', '--sfs_file', default = "None", help = "input site frequency spectrum file")
tiering_group.add_argument('--tiering_only', action='store_true', help='just do tiering (only annotates what is needed for tiering)') # TODO finsh implementing selective annotation

# pgx (pharmgkb/clinvar annotations)
pgx_group = parser.add_argument_group('Pharmacogenomics', 'Arguments for pharmacogenomics module')
pgx_group.add_argument('-n', '--num_threads', default = 1)
pgx_group.add_argument('-b', '--input_bam', default = "None", help = 'input bam file')
pgx_group.add_argument('--reference_sequence', help = 'reference sequence (FASTA file) to accompany input bam')
pgx_group.add_argument('--ref', help='reference (e.g. hg19 or GRCh37)', default='GRCh37')
pgx_group.add_argument('--pgx_only', action='store_true', help='just do pgx (takes care of required annotation, etc.)')
pgx_group.add_argument('--dbsnp', help='path to dbsnp file (used by GATK)')

#debug modes
debug_group = parser.add_argument_group('Debugging parameters')
debug_group.add_argument('--debug_point_annotations', action='store_true', help='enable debug mode for point annotations')
debug_group.add_argument('--point_annotations_only', action='store_true', help='only do point annotations (no functional/range annotations, tiering, etc.)')
debug_group.add_argument('--print_sql_query_only', action='store_true', help='just print out SQL query')
debug_group.add_argument('--skip_join_checks', action='store_true', help='skips join checks for snpeff, annovar, point annotation, range annotation. Without this the log file may get fairly large.')


args = parser.parse_args()

#print help and exit if no arguments are supplied.
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)

# parse YAML
yaml_commands= yaml_utils.parse_yaml_input_files(args.yaml, args.modules)
if(args.download_datasets_only):
	stmp_annotation_util.downloadDBs(yaml_commands, args.dataset_output_dir, args.log)
	sys.exit(0)
elif args.check_datasets_only:
	stmp_annotation_util.checkDBs(yaml_commands, args.dataset_output_dir)
	sys.exit(0)

#open the connection to the database
db_conn = stmp_annotation_util.connect_db(db_file=args.database_file, host_name='', user_name='', db_name='', unix_socket_loc='')

if args.drop_samples:
	stmp_annotation_util.drop_samples(db_conn)

if args.update_db: # check if database setup is requested
	if args.input_directory != None:
		args.update_db = stmp_annotation_util.root_or_cwd(args.input_directory) # complete the filepath if an absolute filepath is not provided.
		stmp_annotation_util.setup_db(args.input_directory, db_conn, not args.force_overwrite) # Launch DB updating process.
	elif args.yaml != None:
		stmp_annotation_util.setup_db_yaml(db_conn, yaml_commands, not args.force_overwrite)
	else:
		print 'Error: neither YAML nor directory with input datasets specified'
		parser.print_help()
		sys.exit(1)

if args.vcf != None: # annotation
	#Files and Directories
	args.vcf = stmp_annotation_util.root_or_cwd(args.vcf) # complete the filepath if an absolute filepath is not provided.
	args.output_dir = stmp_annotation_util.root_or_cwd(args.output_dir) # ditto
	args.scratch_dir = os.path.join(args.output_dir, 'scratch')
	if not os.path.exists(args.scratch_dir):
		os.makedirs(args.scratch_dir)
	
	#Convert multiallelic to single line for easier merging of functional annotations later
	if args.reuse_multiallelic and not args.skip_multiallelic:
		noMultialllic_vcf = stmp_annotation_util.splitMultiallelic(args.vcf, args.scratch_dir, skip_if_exists=True)
		args.vcf = noMultialllic_vcf
	elif not args.skip_multiallelic:
		noMultialllic_vcf = stmp_annotation_util.splitMultiallelic(args.vcf, args.scratch_dir)
		args.vcf = noMultialllic_vcf
	else:
		print 'Skipping multiallelic check'
	
	# strip chr prefix
	stripChr_vcf = stmp_annotation_util.stripChrPrefix(args.vcf, args.scratch_dir, skip_if_exists=False)
	args.vcf=stripChr_vcf

############### MAIN FUNCTIONS (ANNOTATION, TIERING, PGX) ################
	
	# ANNOTATION MAIN FUNCTION
	def annotate(args):		
		if args.region_annotations_only:
			#region annotation does not require sample vcf to be uploaded to our db
			region_outfile = stmp_annotation_util.annotate_range(db_conn, args.vcf, args.scratch_dir, modules_yaml_path=args.modules, yaml_commands=yaml_commands, force_overwrite_beds=args.clean_beds)
			exit(0)
		
		#Functional variant effect prediction
		#snpeff
		if(not args.point_annotations_only):
			[snpeff_proc, sample_name, snpeff_vcf_outfile] = stmp_annotation_util.snpeff(args.vcf, args.scratch_dir) # launch snpEff on the VCF.  This will chug in the background, and the script will wait for the completion of the returned process at the end if it hasn't terminated by that time.
		
		#Actual Annotation
		# upload vcf
		sample_db_path = stmp_annotation_util.upload_vcf(db_conn, args.vcf, args.scratch_dir, args.force_input_overwrite)

		# annovar + region annotation
		if not args.point_annotations_only:
			[annovar_proc, annovar_outfile] = stmp_annotation_util.annovar_annotate_functional(args.vcf, args.scratch_dir)
			region_outfile = stmp_annotation_util.annotate_range(db_conn, args.vcf, args.scratch_dir, modules_yaml_path=args.modules, yaml_commands=yaml_commands, force_overwrite_beds=args.clean_beds) # Find annotations which cover a genomic interval.  This is done with BEDtools.
		
		# point annotation
		point_outfile = stmp_annotation_util.annotate_point(db_conn, args.vcf, args.scratch_dir, sample_db_path, debug=args.debug_point_annotations) # Find annotations which are associated with a single locus.  This is done with a SQL join.
		if args.point_annotations_only:
			exit(0) # stop after point annotation done
	
		#wait for snpeff and annovar to complete, if they haven't already.
		print 'Waiting for snpeff to complete...'
		snpeff_proc.wait()
		print 'Waiting for annovar to complete...'
		annovar_proc.wait()
			
		# convert snpeff vcf to tsv (remove all lines with '#')
		snpeff_tsv_outfile = stmp_annotation_util.snpeff2tsv(sample_name, snpeff_vcf_outfile, args.scratch_dir)
	
		#join the results into a single file
		joined_outfile = stmp_annotation_util.joinFiles(args.vcf, snpeff_tsv_outfile, annovar_outfile, region_outfile, point_outfile, args.output_dir, args.skip_join_checks, yaml_commands=yaml_commands)
		
		# TODO may need to clean up temporary range annotation files (in output/scratch dir) to avoid issues with region annotation.
		# for now (must do this to avoid issues with region annotation)
# 		print 'Cleaning up temporary range annotation files'
# 		cmd = 'rm -f /tmp/stmp2/intersected/*'
# # 		cmd = 'mv /tmp/stmp2/intersected/* /tmp/stmp2/intersected/old/'
# 		subprocess.Popen(cmd, shell=True).wait()
		
		print 'Done annotating ' + joined_outfile
		return joined_outfile
	
	# end annotation function
	
	
	#######################################################################
	####### Variant Tiering ##############
	
	## HELPER FUNCTIONS
	#wrapper to call mapped processes
	def call_process(command):
	    print "Processing command: "+command
	    status = "Return code for command "+command+":"+str(os.system(command))
	    return status
	   
	#function to add rsid to GATK output
	def add_rsid(intervals, in_file, out_file):
		f1 = open(intervals, "r")
		f2 = stmp_annotation_util.open_compressed_or_regular(in_file, "r")
		f3 = open(out_file, "w")
	
		rsdict = {}
		for line in f1:
	        	if ("#" in line) == 0:
	                	linelist = line.split("\t")
	                	rsdict[linelist[0]+":"+linelist[1]] = linelist[2].replace("target_", "").replace("\n", "")
	
		while 1:
			temp = f2.readline()
	    		if not temp:
	       			break
	    		else:
	        		if (("#" in temp) == 0):
	                		linelist = temp.split("\t")
	                		if rsdict.has_key(linelist[0]+":"+linelist[1]):
	                        		f3.write(linelist[0]+"\t"+linelist[1]+"\t"+rsdict[linelist[0]+":"+linelist[1]]+"\t"+"\t".join(linelist[3:len(linelist)]))
	        		else:
	                		f3.write(temp)
		f1.close()
		f2.close()
		f3.close()
	
	
	# TIERING (VARIANT PRIORITIZATION) MAIN FUNCTION
	def tier(args, annotation_joined_outfile, yaml_cmds, output_dir=None, tier_name='Global'):
		## MAIN CODE
		print 'Performing variant prioritization'
		
		if(output_dir == None):
			output_dir = os.path.join(args.output_dir, 'tiering_allvars')
		else: 
			output_dir = os.path.join(output_dir, 'tiering_allvars')
			
		# targeted tiering only if BAM file provided? (TODO)
		# stmp_tiering_util.tiers_target(os.path.join(args.output_dir, 'tiers_target.tsv'), os.path.join(args.output_dir, 'tiers_target'), args.target_genes, pop=args.ethnicity)
		
		# standard VCF prioritization (tiering)
		stmp_tiering_util.tiers_allvars(annotation_joined_outfile, output_dir, args.target_genes, pop=args.ethnicity, yaml_cmds=yaml_cmds)
		
		# SFS filtering (TODO)
# 		if args.sfs_file != "None":	
# 			for i in range(1,5):
# 				stmp_tiering_util.filter_sfs(os.path.join(args.output_dir, "allvars.tier"+str(i)+".txt"), args.sfs_file, os.path.join(args.output_dir, "allvars.tier"+str(i)+"-sfs_filtered.txt"), 2)
	
	# end tiering function
	
	# real tiering function (uses above function as helper)
	def tier_real(args, joined_outfile, yaml_commands):
		#tiering is separate
		tiering_output_dirs = []
		
		# 1. candidate genes (user-specified target gene list)
		if(args.target_genes != None):
			args.target_genes = general_utils.root_or_cwd(args.target_genes)
			candidate_out_dir = os.path.join(args.output_dir, 'Candidate')
			if(not os.path.isdir(candidate_out_dir)):
				os.makedirs(candidate_out_dir)
			tier(args, joined_outfile, yaml_cmds=yaml_commands, output_dir=candidate_out_dir)
			tiering_output_dirs.append(candidate_out_dir)
		
		# 2. global
		global_tiering_out_dir = os.path.join(args.output_dir, 'Global')
		if(not os.path.isdir(global_tiering_out_dir)):
			os.makedirs(global_tiering_out_dir)
		args.target_genes = None # forces global tiering (no filtering based on a certain gene list)
		tier(args, joined_outfile, yaml_cmds=yaml_commands, output_dir=global_tiering_out_dir)
		tiering_output_dirs.append(global_tiering_out_dir)
		
		#3. other gene lists (specified in YAML)
		tiering_gene_lists = yaml_commands[yaml_keys.kModules][yaml_keys.kTiering][yaml_keys.kTTargetGeneLists]
		for tiering_gene_list in tiering_gene_lists:
			out_dir = os.path.join(args.output_dir, tiering_gene_list)
			if(not os.path.isdir(out_dir)):
				os.makedirs(out_dir)
			args.target_genes = yaml_utils.get_abs_path(tiering_gene_lists[tiering_gene_list])
			tier(args, joined_outfile, yaml_cmds=yaml_commands, output_dir=out_dir)
			tiering_output_dirs.append(out_dir)
		
		
		# Generate final output as excel workbook
		final_tiering_out_file_path = stmp_tiering_util.tiers2xls(tiering_output_dirs, args.output_dir)
		print '**** Tiering output written to ' + str(final_tiering_out_file_path) + ' *****'
		return final_tiering_out_file_path
	
	#################################################################################################
	
	########### Pharmgkb and ClinVar Annotation #######################

	# PGX MAIN FUNCION (PHARMGKB AND CLINVAR ANNOTATION)
	def pgx(args, output_dir=None):
		# consts
		if(output_dir == None):
			args.output = os.path.join(args.output_dir, 'pgx')
		else:
			args.output = output_dir
		
		if(not os.path.isdir(args.output)):
			os.makedirs(args.output)
		
		########## MAIN CODE #####################
	
		print str(datetime.datetime.now()) + ': Performing PharmGKB + ClinVar (pgx) annotation'

		#pharmgkb annotation
		pgxUtils.pgx_annotator(args.vcf, os.path.join(resources_path, "pgx_vars", "clinical_ann_metadata-snvs.txt"), os.path.join(args.output, "pharmacogenomics"))
		pgxUtils.star_caller(os.path.join(resources_path, "pgx_haps/"), args.vcf, os.path.join(args.output, "pharmacogenomics"))

		print str(datetime.datetime.now()) + ': Done with pgx/clinvar annotation'
		
	# end pgx function
	
	# to be tested
	def write_output_yaml(yaml_commands, output_dir):
		outfile = open(os.path.join(output_dir, 'config.yml'), 'w')
		yaml.dump(yaml_commands, outfile, default_flow_style=False)
		return outfile

###################### MAIN CODE (logic to call main functions) ###############

	#just print sql query if specified
	if args.print_sql_query_only:
		sample_db_path = stmp_annotation_util.get_sample_db_path(args.scratch_dir, stmp_annotation_util.getSampleName(args.vcf))
		stmp_annotation_util.annotate_point(db_conn, args.vcf, args.scratch_dir, sample_db_path, print_sql_query_only=True, debug=args.debug_point_annotations)
		exit(0)
	
	if args.print_region_cmd_only:
		stmp_annotation_util.annotate_range(db_conn, args.vcf, args.scratch_dir, modules_yaml_path=args.modules, yaml_commands=yaml_commands, print_range_cmd_only=True, force_overwrite_beds=args.clean_beds)
		exit(0)
	
	# pgx doesn't depend on other annotations
	if(args.pgx_only):
		pgx(args)
		exit(0)
	
	# tiering requires annotation to already be done
	if(args.tiering_only):
		joined_outfile = stmp_annotation_util.generateJoinedOutfilePath(args.output_dir, stmp_annotation_util.getSampleName(args.vcf)) # generates path to outfile but does not create it
		stmp_annotation_checker.check_annotated_output(joined_outfile)
		tier_real(args, joined_outfile, yaml_commands)
		exit(0)
	
	if(args.annotate_only):
		joined_outfile = annotate(args)
		stmp_annotation_checker.check_annotated_output(joined_outfile)
		exit(0)
		
		
	######## DEFAULT #######
	
	# Run full pipeline (global, clinical, candidate, secondary, pgx)
	
	#annotation is common to all modes
	joined_outfile = annotate(args)
	#check annotated output file
	stmp_annotation_checker.check_annotated_output(joined_outfile)
	
	# tier
	tier_real(args, joined_outfile, yaml_commands)
	
	# 5. pgx - goes in its own folder by default
	pgx(args)
	
	
	#remove sample from database if desired
	if(args.drop_sample):
		stmp_annotation_util.drop_sample(db_conn, args.vcf)
	
	# clean up samples in database if needed
	if(args.drop_samples):
		stmp_annotation_util.drop_samples(db_conn)
	
