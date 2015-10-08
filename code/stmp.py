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

stmp.py --update_db --yaml=config/annotation.yml # loads all datasets specified in the YAML file into our database

stmp.py --vcf={input vcf or vcf.gz file} --output_dir={output directory} # applies annotation and tiering to the input VCF, outputting it in the specified output directory.

For a complete list of parameters, just run stmp.py with no options.


'''


import argparse
import stmp_annotation_util
import stmp_tiering_util
import pgxUtils
import multiprocessing
from multiprocessing import pool
import time
from time import strftime, localtime
import os
import sys
import yaml
import subprocess
import datetime
import gzip

cur_path = os.path.dirname(os.path.abspath(__file__))
resources_path= os.path.join(cur_path, os.path.pardir, 'data')

# parse commandline options.
parser = argparse.ArgumentParser(description = 'Sequence to Medical Phenotypes')

parser.add_argument("--config", dest='yaml', help = 'yaml-formatted annotation/db update specification', default=os.path.join(cur_path, 'config', 'datasets.yml'))
# parser.add_argument('--modules', help='yaml-formatted options for each module (annotation, tiering, pgx, trio)', default=os.path.join(cur_path, 'config', 'modules.yml'))

# arguments related to updating the database
db_update_group = parser.add_argument_group("Database initialization/update", 'Initialize DB, or update if already extant.')
db_update_group.add_argument("--update_db", action='store_true', 
	help = 'Initialize or update MySQL database from text annotation files.  Please provide either YAML configuration file or directory containing annotation database text files (see below parameters)')
db_update_group.add_argument('--input_directory', help='Directory containing annotation database text files.')
db_update_group.add_argument("--skip", help = 'Skip data sources already present in the database', default = False, const = True, nargs = '?')
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
annotate_group.add_argument('--db_bed_dir', help='location of BED files generated for range annotation (default: ../db/db_beds)', default=os.path.join(cur_path, os.pardir, 'db', 'db_beds/'))
annotate_group.add_argument('--annotate_only', action='store_true', help='just do annotation (no tiering, pgx, etc.)')
annotate_group.add_argument('--drop_samples', action='store_true', help='drops sample tables (use this before distributing the database)')
annotate_group.add_argument('--drop_sample', action='store_true', help='drops this sample from the database after completing annotation')
annotate_group.add_argument('--region_annotations_only', action='store_true', help='just do region annotation')
annotate_group.add_argument('--print_region_cmd_only', action='store_true', help="just print region annotation command (don't actually do the annotation or anything else")
# annotate_group.add_argument('--sorted_vcf', action='store_true', help='use this if the input VCF is presorted by chr, then start position. Range annotation will then use more memory-efficient algorithms.')

# variant tiering (prioritization)
tiering_group = parser.add_argument_group('Variant Prioritization', 'Arguments for variant prioritization')
tiering_group.add_argument('-t', '--target_genes', default = os.path.join(resources_path, 'clinvar', 'clin_var_curated_genes.txt'),
help = 'target gene file, default clinvar') # 						cur_path+"/resource/clinvar/clin_var_curated_genes.txt", help = 'target gene file, default clinvar')
tiering_group.add_argument('-e', '--ethnicity', default = "CEU", help = 'ethnicity for allele frequency filtering')
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
# debug_group.add_argument('--compact_logs')
debug_group.add_argument('--skip_join_checks', action='store_true', help='skips join checks for snpeff, annovar, point annotation, range annotation. Without this the log file may get fairly large.')


args = parser.parse_args()

#print help and exit if no arguments are supplied.
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)

# parse YAML if specified
if(args.yaml != None):
	yaml_commands = stmp_annotation_util.parse_yaml(args.yaml)
	if(args.download_datasets_only):
		stmp_annotation_util.downloadDBs(yaml_commands, args.dataset_output_dir, args.log)
		sys.exit(0)
	elif args.check_datasets_only:
		stmp_annotation_util.checkDBs(yaml_commands, args.dataset_output_dir)
		sys.exit(0)
	

# clean up any existing bed files if needed - NOTE: no longer needed
# if args.clean_beds:
# 	print 'Cleaning up old BED files'
# 	cmd = 'rm -f /tmp/stmp2/db_beds/*'
# 	#cmd = 'mv /tmp/stmp2/db_beds/* /tmp/stmp2/db_beds/old/'
# 	subprocess.Popen(cmd, shell=True).wait()

#open the connection to the database
# if(args.database_file != None):
db_conn = stmp_annotation_util.connect_db(db_file=args.database_file, host_name='', user_name='', db_name='', unix_socket_loc='')
# else:
# 	db_conn = stmp_annotation_util.connect_db(host_name='', user_name='', db_name='', unix_socket_loc='')

if args.drop_samples:
	stmp_annotation_util.drop_samples(db_conn)

if args.update_db: # check if database setup is requested
	if args.input_directory != None:
		args.update_db = stmp_annotation_util.root_or_cwd(args.input_directory) # complete the filepath if an absolute filepath is not provided.
		stmp_annotation_util.setup_db(args.input_directory, db_conn, args.skip) # Launch DB updating process.
	elif args.yaml != None:
		stmp_annotation_util.setup_db_yaml(db_conn, yaml_commands, args.skip)
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
	
	########################################################################################################################
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
			region_outfile = stmp_annotation_util.annotate_range(db_conn, args.vcf, args.scratch_dir, db_bed_dir=args.db_bed_dir)
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
			region_outfile = stmp_annotation_util.annotate_range(db_conn, args.vcf, args.scratch_dir, db_bed_dir=args.db_bed_dir) # Find annotations which cover a genomic interval.  This is done with BEDtools.
		
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
		joined_outfile = stmp_annotation_util.joinFiles(args.vcf, snpeff_tsv_outfile, annovar_outfile, region_outfile, point_outfile, args.output_dir, args.skip_join_checks)

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
	def tier(args, annotation_joined_outfile, yaml_cmds):
		## MAIN CODE
		print 'Performing variant prioritization'
		
		# targeted tiering only if BAM file provided?
		# stmp_tiering_util.tiers_target(os.path.join(args.output_dir, 'tiers_target.tsv'), os.path.join(args.output_dir, 'tiers_target'), args.target_genes, pop=args.ethnicity)
		
		# standard VCF prioritization (tiering)
		stmp_tiering_util.tiers_allvars(joined_outfile, os.path.join(args.output_dir, 'tiering_allvars'), args.target_genes, pop=args.ethnicity, yaml_cmds=yaml_cmds)
		
		if args.sfs_file != "None":	
			for i in range(1,5):
				stmp_tiering_util.filter_sfs(os.path.join(args.output_dir, "allvars.tier"+str(i)+".txt"), args.sfs_file, os.path.join(args.output_dir, "allvars.tier"+str(i)+"-sfs_filtered.txt"), 2)
		
		
		print("Done with variant prioritization. Check output files in " + args.output_dir)
	
	# end tiering function
	
	#################################################################################################
	
	########### Pharmgkb and ClinVar Annotation #######################

	# PGX MAIN FUNCION (PHARMGKB AND CLINVAR ANNOTATION)
	def pgx(args):
		# consts
		args.output = args.output_dir
		
		########## HELPER METHODS ###############
		
		#wrapper to annotate snvs/indels
	# 	def annotate_snvs(snv_in, snv_out_stem):
	# 	
	# 		#convert input file to annovar format
	# 		if ".vcf" in snv_in:
	# 			os.system("perl "+cur_path+"/stanovar/convert2annovar.pl "+snv_in+" --includeinfo -format vcf4 > "+snv_out_stem+".snv.annovar.txt")
	# 		else:
	# 			print >> sys.stderr, "Error in annotate_vars.py: Invalid vcf input file format"
	# 			exit(1)
	# 	
	# 		#read in command file
	# 		#yaml_in = "/stanovar/config/annotation_testing_subset.yaml"
	# 		#print("USING TESTING ANNOTATION SUBSET")
	# 		yaml_in = "/stanovar/config/annotation.yaml"
	# 		commands = parse_yaml_commands(snv_out_stem, yaml_in)
	# 		return commands
		
		########## MAIN CODE #####################
	
		print str(datetime.datetime.now()) + ': Performing PharmGKB + ClinVar (pgx) annotation'
		
# 		# GENOTYPING
# 		#pool of processors for concurrent processing
# 		pool = multiprocessing.Pool(processes = int(args.num_threads))
# 
# 		command_call = []
# 		#call targeted genotypes and perform clinvar and PharmGKB annotation, if BAM file provided
# 		if args.input_bam != "None":
# 			if ".bam" in args.input_bam:
# 				# choose appropriate resource files based on specified reference (hg19, GRCh37, etc.). Default file (e.g. pgx_cd_5.7.13.snvs.vcf.gz) is grch37 (does not include chr prefix), other file (pgx_cd_5.7.13_hg19.snvs.vcf.gz) is hg19 (includes chr prefix).
# 				pgx_snv_interval_filename = 'pgx_cd_5.7.13.snvs.vcf.gz_b37.vcf.gz' if args.ref.lower() == 'grch37' else 'pgx_cd_5.7.13.snvs_hg19.vcf.gz' if args.ref.lower() == 'hg19' else ''
# 				pgx_indel_interval_filename = 'pgx_cd_5.7.13.b37.indels.vcf' if args.ref.lower() == 'grch37' else 'pgx_cd_5.7.13.indels.vcf' if args.ref.lower() == 'hg19' else ''
# # 				clinvar_indel_interval_filename = 'clinvar_3.15.13.b37.indels.vcf' if args.ref.lower() == 'grch37' else 'clinvar_3.15.13.indels.vcf' if args.ref.lower() == 'hg19' else ''
# 				#debug
# # 				clinvar_indel_interval_filename = 'clinvar_3.15.13.indels.vcf'
# 				clinvar_indel_interval_filename = 'clinvar_indels.vcf.recode.vcf'
# 				
# 				# then call appropriate commands
# 				cmd = "sh "+cur_path+"/target_snvs.sh "+args.input_bam+" "+cur_path+"/resource/intervals/clinvar.interval_list "+args.output+"/clinvar.snvs"+' '+args.reference_sequence+' '+args.dbsnp
# 				#debug
# 				print 'cmd: ' + cmd
# 				command_call.append(cmd)
# 				
# 				cmd = "sh "+cur_path+"/target_snvs.sh "+args.input_bam+" "+cur_path+"/resource/intervals/"+pgx_snv_interval_filename+" "+args.output+"/pgx.snvs"+' '+args.reference_sequence+' '+args.dbsnp
# 				#debug
# 				print 'cmd: ' + cmd
# 				command_call.append(cmd)
# 				
# 				cmd = "sh "+cur_path+"/target_indels.sh "+args.input_bam+" "+cur_path+"/resource/intervals/"+clinvar_indel_interval_filename+" "+args.output+"/clinvar.indels"+' '+args.reference_sequence+' '+args.dbsnp
# 				#debug
# 				print 'cmd: ' + cmd
# 				command_call.append(cmd)
# 				
# 				cmd = "sh "+cur_path+"/target_indels.sh "+args.input_bam+" "+cur_path+"/resource/intervals/"+pgx_indel_interval_filename+" "+args.output+"/pgx.indels"+' '+args.reference_sequence+' '+args.dbsnp
# 				#debug
# 				print 'cmd: ' + cmd
# 				command_call.append(cmd)
# 
# 				#log file specification
# 				log_file = open(args.output+"/stmp2.log", "w")
# 				log_file.write(">>> Sequence to medical phenotypes log file <<<\n")
# 				log_file.write("\nstmp2 started: "+strftime("%a, %d %b %Y %H:%M:%S +0000", localtime())+"\n\n")
# 
# 				#map calling commands
# 				res_calls = pool.map(call_process, command_call)
# 				for item in res_calls:
# 					log_file.write(item+"\n")
# 
# 				#merge clinvar, and pgx calls and re-annotate rsid
# 				add_rsid(cur_path+"/resource/intervals/"+pgx_snv_interval_filename, args.output+"/pgx.snvs.filtered.vcf", args.output+"/pgx.snvs.annotated.vcf")
# 				add_rsid(cur_path+"/resource/intervals/"+pgx_indel_interval_filename, args.output+"/pgx.indels.filtered.vcf", args.output+"/pgx.indels.annotated.vcf")
# 	#			add_rsid(cur_path+"/resource/intervals/clinvar_3.15.13.snvs.vcf", args.output+"/clinvar.snvs.filtered.vcf", args.output+"/clinvar.snvs.annotated.vcf")
# 	#			add_rsid(cur_path+"/resource/intervals/clinvar_3.15.13.indels.vcf", args.output+"/clinvar.indels.filtered.vcf", args.output+"/clinvar.indels.annotated.vcf")
# 				os.system("sh "+cur_path+"/combine_variants.sh "+args.output+"/pgx.snvs.annotated.vcf "+args.output+"/pgx.indels.annotated.vcf "+args.output+"/pgx.all.vcf")
# 				os.system("sh "+cur_path+"/combine_variants.sh "+args.output+"/clinvar.snvs.filtered.vcf "+args.output+"/clinvar.indels.filtered.vcf "+args.output+"/clinvar.all.vcf")

				#annotate and prioritize clinvar variants using Stanovar
	# 			command_ann_clin = annotate_snvs(args.output+"/clinvar.all.vcf", args.output+"/clinvar")
	# 			res_ann = pool.map(call_process, command_ann_clin)
	# 			for item in res_ann:
	#         			log_file.write(item+"\n")
	# 			os.system("perl "+cur_path+"/stanovar/summarize_annovar.pl "+args.output+"/clinvar"+" "+cur_path+"/stanovar/humandb '"+get_head(args.output+"/clinvar.all.vcf")+"'")
	# 			diseaseUtils.tiers_target(args.output+"/clinvar.genome_summary.tsv", args.output+"/clinvar", args.target_genes, 0.01, args.ethnicity)
	# 			if args.sfs_file != "None":	
	# 				for i in range(1,5):
	# 					diseaseUtils.filter_sfs(args.output+"/clinvar.tier"+str(i)+".txt", args.sfs_file, args.output+"/clinvar.tier"+str(i)+"-sfs_filtered.txt", 2)

		#pharmgkb annotation
		
		pgxUtils.pgx_annotator(args.vcf, os.path.join(resources_path, "pgx_vars", "clinical_ann_metadata-snvs.txt"), os.path.join(args.output, "pharmacogenomics"))
		pgxUtils.star_caller(os.path.join(resources_path, "pgx_haps/"), args.vcf, os.path.join(args.output, "pharmacogenomics"))
		
# 		pgxUtils.pgx_annotator(os.path.join(args.output, "pgx.all.vcf"), os.path.join(cur_path, "resource/pgx_vars/clinical_ann_metadata-snvs.txt"), os.path.join(args.output, "pharmacogenomics"))
# 		pgxUtils.star_caller(os.path.join(cur_path, "resource/pgx_haps/"), os.path.join(args.output, "pgx.all.vcf"), os.path.join(args.output, "pharmacogenomics"))			
# 	else:
# 		print >> sys.stderr, "Error in stmp2.py - input file does not appear to be bam format"
# 		exit(1)

		print str(datetime.datetime.now()) + ': Done with pgx/clinvar annotation'
		
	# end pgx function


###################### MAIN CODE (logic to call main functions) ###############

	#just print sql query if specified
	if args.print_sql_query_only:
		sample_db_path = stmp_annotation_util.get_sample_db_path(args.scratch_dir, stmp_annotation_util.getSampleName(args.vcf))
		stmp_annotation_util.annotate_point(db_conn, args.vcf, args.scratch_dir, sample_db_path, print_sql_query_only=True, debug=args.debug_point_annotations)
		exit(0)
	
	if args.print_region_cmd_only:
		stmp_annotation_util.annotate_range(db_conn, args.vcf, args.scratch_dir, db_bed_dir=args.db_bed_dir, print_range_cmd_only=True)
		exit(0)
	
	# pgx doesn't depend on other annotations
	if(args.pgx_only):
		pgx(args)
		exit(0)
	
	# tiering requires annotation to already be done
	if(args.tiering_only):
		joined_outfile = stmp_annotation_util.generateJoinedOutfilePath(args.output_dir, stmp_annotation_util.getSampleName(args.vcf)) # generates path to outfile but does not create it
		tier(args, joined_outfile, yaml_cmds=yaml_commands)
		exit(0)
	
	if(args.annotate_only):
		annotate(args)
		exit(0)
	
	# by default, run full pipeline
	joined_outfile = annotate(args)
	tier(args, joined_outfile, yaml_cmds=yaml_commands)
	
	# comment out the below line for stable branch (currently doesn't work)
	pgx(args)
	
	
	#remove sample from database if desired
	if(args.drop_sample):
		stmp_annotation_util.drop_sample(db_conn, args.vcf)
	
	# clean up samples in database if needed
	if(args.drop_samples):
		stmp_annotation_util.drop_samples(db_conn)
	
