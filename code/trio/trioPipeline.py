#Rick Dewey 6.29.12
#Modified by James Priest 8.26.13
#Python code for Trio case where offspring is affected
#usage: trioPipeline.py <input vcf> <output file stem> <annovar> <matrix_path> <offspringID> <fatherID> <motherID> <optional commands...>
#!/usr/bin/env python

import os
import re
import fileinput
import sys
import getopt
import argparse
import math
import vcfUtils
import hmmUtils
import numpy as np
import math

  
#get args from command line
parser = argparse.ArgumentParser(description = 'Finds variants consistent with sporadic disease or autosomal recessive transmission in trios')
parser.add_argument('input', help = 'input vcf file')
parser.add_argument('output', help = 'file stem for formatted output')
parser.add_argument('annovar', help = 'annovar environment')
parser.add_argument('matrix_path', help = 'path to emission and transition matrices')
parser.add_argument('offspringID', help = 'ID for offspring')
parser.add_argument('fatherID', help = 'ID for father')
parser.add_argument('motherID', help = 'ID for mother')
parser.add_argument('-v','--annotate_variants', help = 'perform annovar annotation, default = TRUE', default = 'TRUE')
parser.add_argument('-r','--reference', help = 'reference coordinates for annovar, default hg19', default = "hg19")
parser.add_argument('-p','--dbsnp', help = 'version of dbsnp, default 135', default = "135")
parser.add_argument('-g','--gene_model', help = 'gene model, default UCSC knowngene', default = "knowngene")
parser.add_argument('-d','--depth', help = 'coverage depth threshold for inclusion, default 0', default = 0)
parser.add_argument('-a','--alt_depth', help = 'threshold for reads containing alt allele, default 0', default = 0)
parser.add_argument('-f','--allele_frequency', help = 'allele frequency cut-off for inclusion, default 0.01', default = 0.01)
parser.add_argument('-b','--background_population', nargs = "+", help = 'background population for allele frequency estimates', default = '1000g2012apr_ALL ESP6500_ALL')
parser.add_argument('-q','--gt_quality', help = 'genotype quality cut-off for inclusion, default 0', default = 0)
parser.add_argument('-m','--mapping_quality', help = 'mapping quality cut-off for inclusion, default 0', default = 0)
parser.add_argument('-z','--mq0', help = 'mapping quality 0 cut-off for inclusion, default 0', default = 0)
parser.add_argument('-l','--AVR_cutoff', help = 'user defined AVR cutoff for inclusion default 0', default = 0.0)
parser.add_argument('-e','--error_prob', help = 'per-base error estimate, default 0.005', default = 0.005)
parser.add_argument('-c','--compression_prob', help = 'emission probability for uniform heterozygosity for compressions, default 0.66', default = 0.66)
parser.add_argument('-i','--mie_prob', help = 'emission probability for MIE in MIE rich regions, default 0.33', default = 0.33)
parser.add_argument('-s','--self_transition_prob', help = 'self-self transition probability, default 0.999999', default = 0.999999)
parser.add_argument('-x','--exome', help = 'data is exome data, default FALSE', default = 'FALSE')
#parser.add_argument('-t','--platform', help = 'sequencing platform ILLUMINA, CG, or RTG', default = 'ILLUMINA')
args = parser.parse_args()

vcf_in = args.input
final_out = args.output
depth = int(args.depth)
altDepth = int(args.alt_depth)
freq = float(args.allele_frequency)
gtQual = float(args.gt_quality)
mapQual = float(args.mapping_quality)
mapQual0 = float(args.mq0)
pl = float(args.AVR_cutoff)
mat_path = args.matrix_path
error = args.error_prob
compression = args.compression_prob
mie = args.mie_prob
selfprob = args.self_transition_prob
annovar = args.annovar
ref = args.reference
dbSNP = args.dbsnp
geneModel = args.gene_model
exome = args.exome
annotate = args.annotate_variants
backpop = args.background_population
platform = "RTG"


#create matrices for HMM, all probabilities in log10 space to avoid underflow
transmat, emitmat, startp, states, inmap = hmmUtils.make_hmm(mat_path, selfprob, compression, mie, error, "t")

#will run annovar annotations if not done already, needs error message for non-vcf file input
if annotate == "TRUE":
    #format vcf file for annovar (all chromosomes)
    print "\nStep 1: Formatting vcf file for annovar\n"
    print "Converting vcf to annovar format...\n"
    os.system("perl "+annovar+"convert2annovar.pl "+vcf_in+" --includeinfo -format vcf4 > "+vcf_in.replace(".vcf", "_annovar.txt"))

    #annotate using modified annovar, annotates all variants at once to avoid slow-downs associated with loading annotation files into memory repeatedly
    print "\nStep 2: Annotating variants using annovar\n"
    os.system("perl "+annovar+"summarize_annovarRDv2.pl "+vcf_in.replace(".vcf", "_annovar.txt ")+annovar+"humandb/ --buildver "+ref+" --verdbsnp "+dbSNP+" --genetype "+geneModel+" --remove")

    #reformat header for hmm and sporadic trio allele finding
    print "\nStep 3: Adding vcf header to annotated vcf file\n"
    temp = open(vcf_in.replace(".vcf", "_annovar.txt.genome_summary.txt"), "r")
    head1 = temp.readline()
    headvcf = vcfUtils.get_head(vcf_in)
    newhead = head1.replace("Otherinfo\n", headvcf)

    annovar_out = open(vcf_in.replace(".vcf", "_genome_summary_formattedALL.txt"), "w")
    annovar_out.write(newhead)
    for line in temp:
        annovar_out.write(line)
    annovar_out.close()
    temp.close()
    os.system("rm "+vcf_in.replace(".vcf", "_annovar.txt.genome_summary.txt"))

#otherwise skips annovar annotation
else:
    print "\nVariants are already annotated. Skipping steps 1-3."

#only run QA hmm on genome data
if exome == "FALSE":
    print "\nPreparing data for QA HMM...\n"
    
    #split annotated vcf file for HMM
    if annotate == "TRUE":
        vcfUtils.splitFiles(vcf_in.replace(".vcf", "_genome_summary_formattedALL.txt"), vcf_in.replace(".vcf", "_genome_summary_formatted"))
    else:
        vcfUtils.splitFiles(vcf_in, vcf_in.replace(".txt", "_genome_summary_formatted"))
    
    #QA hmm
    print "\nStep 4: Running QA HMM"

    chrom_arr = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]

    for chromh in chrom_arr:
        print "\nWorking on chromosome: "+str(chromh)
        if annotate == "TRUE":
            fin = open(vcf_in.replace(".vcf", "_genome_summary_formatted")+"_chr"+chromh+".txt", "r")
        else:
            fin = open(vcf_in.replace(".txt", "_genome_summary_formatted")+"_chr"+chromh+".txt", "r")
    
        header = fin.readline()
        #grab proper columns from header for variants, need error message for improper input
        headlist = header.replace("\n", "").split("\t")
        offspring = headlist.index(args.offspringID)
        father = headlist.index(args.fatherID)
        mother = headlist.index(args.motherID)
    
        obs_seq = []
        position = []
        res_dict = {}

        counterDN = 0
        counterHZGC = 0
        counter = 0

        #counter for number of informative allele assortments for determing paternity/maternity
        inform_parents = 0
        non_paternal = 0
        non_maternal = 0

        print "Checking sample consistency and creating observation array..."

        for line in fin.readlines():
            counter+=1
            linelist = line.split("\t")
            chrom = linelist[headlist.index("CHROM")]
            posit = linelist[headlist.index("POS")]
            info = linelist[headlist.index("INFO")]
            #format = linelist[headlist.index("FORMAT")]
            mq, mq0 = vcfUtils.parse_info(info)
    
            if ("#" in line) == 0 and mq >= mapQual and mq0 <= mapQual0 and ("PASS" in line):
                #stat_locations = vcfUtils.stat_grabber(linelist
                geno_father = vcfUtils.allele_coder(linelist[father], depth, altDepth, gtQual, pl, platform)
                geno_mother = vcfUtils.allele_coder(linelist[mother], depth, altDepth, gtQual, pl, platform)
                geno_offspring = vcfUtils.allele_coder(linelist[offspring], depth, altDepth, gtQual, pl, platform)
        
                alleles2 = geno_father+","+geno_mother+","+geno_offspring
                if ((("3" in alleles2) == 0) and (("4" in alleles2) == 0)):
                    counter+=1

                    #checks for paternity, maternity
                    if geno_father == "0" and geno_mother == "2":
                        inform_parents += 1
                        if geno_offspring == "2":
                            non_paternal += 1
                        if geno_offspring == "0":
                            non_maternal += 1

                    elif geno_father == "2" and geno_mother == "0":
                        inform_parents += 1
                        if geno_offspring == "0":
                            non_paternal += 1
                        if geno_offspring == "2":
                            non_maternal += 1
                    
                    #simple de novo case
                    if alleles2 == "0,0,1":
                        counterDN+=1

                    #hemizygosity or gene conversion case
                    elif alleles2 == "0,0,2" or alleles2 == "1,0,2" or alleles2 == "0,1,2" or alleles2 == "1,2,0" or alleles2 == "2,1,0":
                        counterHZGC+=1

                    #creates an array of clean observations for HMM
                    if inmap.has_key(alleles2):
                        obs_seq.append(int(inmap[alleles2]))
                        position.append(posit)

        fin.close()
                
        cent_np = float(non_paternal)/float(inform_parents)
        cent_nm = float(non_maternal)/float(inform_parents)
        cent_DN = float(counterDN+counterHZGC)/float(counter)

        print "...Done"

        print "\nNon-transmitted paternal alleles (%): "+str(cent_np*float(100))
        print "Non-transmitted maternal alleles (%): "+str(cent_nm*float(100))
        print "Novel alleles (%): "+str(cent_DN*float(100))

        #warning for high rates of non-transmission
        if cent_np > 0.05:
            print "\nWARNING! High rate of non-transmission of paternal alleles: consider sample mix-up or non-paternity...\n"
        if cent_nm > 0.05:
            print "\nWARNING! High rate of non-transmission of maternal alleles: consider sample mix-up or non-maternity...\n"
        if cent_DN > 0.05:
            print "\nWARNING! High rate of novel alleles in child: consider sample mix-up or non-maternity or non-paternity...\n"

        print "\nRunning HMM..."

        #run Viterbi algorithm on HMM defined above
        vit = hmmUtils.viterbi(np.array(obs_seq), states, startp, transmat, emitmat)
        vit_path = vit[1]

        #bind viterbi path to position dict for annotation
        for i in range(0, len(vit_path)):
            res_dict[position[i]] = vit_path[i]

        print "...Done"
        
        if annotate == "TRUE":
            print "\nWriting to output file "+vcf_in.replace(".vcf", "_chr"+chromh+"_annotated_hmmfiltered.txt")
            fin = open(vcf_in.replace(".vcf", "_genome_summary_formatted")+"_chr"+chromh+".txt", "r")
            fout = open(vcf_in.replace(".vcf", "_chr"+chromh+"_annotated_hmmfiltered.txt"), "w")
        else:
            print "\nWriting to output file "+vcf_in.replace(".txt", "_chr"+chromh+"_annotated_hmmfiltered.txt")
            fin = open(vcf_in.replace(".txt", "_genome_summary_formatted")+"_chr"+chromh+".txt", "r")
            fout = open(vcf_in.replace(".txt", "_chr"+chromh+"_annotated_hmmfiltered.txt"), "w")
            
        header = fin.readline()
        fout.write(header)

        #temp variables to allow bridging alleles not used in HMM
        firstpath = vit_path[0]
        lastpath = firstpath

        for line in fin.readlines():
    
            if ("#" in line) == 0:
                linelist = line.split("\t")
                chrom = linelist[headlist.index("CHROM")]
                posit = linelist[headlist.index("POS")]
                info = linelist[headlist.index("INFO")]
                format = linelist[headlist.index("FORMAT")]
                filter_field = linelist[headlist.index("FILTER")]

                #write filter flag for MIE/compression if it is present, none if otherwise
                if res_dict.has_key(posit):
                    if res_dict[posit] == 0:
                        fout.write(line.replace(filter_field, "MIEregion"))
                    elif res_dict[posit] == 1:
                        fout.write(line.replace(filter_field, "CompressionRegion"))
                    elif res_dict[posit] == 2:
                        fout.write(line)
                    else:
                        print res_dict[posit]
                    lastpath = res_dict[posit]

                #uses filter from last position if not used in HMM
                else:
                    if lastpath == 0:
                        if ("PASS" in line) == 0:
                            fout.write(line.replace(filter_field, filter_field+";MIEregion"))
                        else:
                            fout.write(line.replace(filter_field, "MIEregion"))
                    elif lastpath == 1:
                        if ("PASS" in line) == 0:
                            fout.write(line.replace(filter_field, filter_field+";CompressionRegion"))
                        else:
                            fout.write(line.replace(filter_field, "CompressionRegion"))
                    elif lastpath == 2:
                        fout.write(line)
                
            else:
                fout.write(line)
        fin.close()

        if annotate == "TRUE":
            os.system("rm "+vcf_in.replace(".vcf", "_genome_summary_formatted")+"_chr"+chromh+".txt")
        else:
            os.system("rm "+vcf_in.replace(".txt", "_genome_summary_formatted")+"_chr"+chromh+".txt")
                  
        fout.close()

    print "\nMerging files...\n"

    #merge annotated file for downstream tools
    if annotate == "TRUE":
        vcfUtils.mergeFiles(vcf_in.replace(".vcf", ""), vcf_in.replace(".vcf", "_annotated_hmmfilteredALL.txt"))
    else:
        vcfUtils.mergeFiles(vcf_in.replace(".txt", ""), vcf_in.replace(".txt", "_annotated_hmmfilteredALL.txt"))

else:
    print "\nExome data: Skipping step 4 (hmm QA)" 

#specific to sporadic trio allele finding
print "\nStep 5: finding alleles consistent with sporadic disease\n"

if exome == "FALSE" and annotate == "TRUE":
    fin = open(vcf_in.replace(".vcf", "_annotated_hmmfilteredALL.txt"), "r")
elif exome == "FALSE" and annotate == "FALSE":
    fin = open(vcf_in.replace(".txt", "_annotated_hmmfilteredALL.txt"), "r")
elif exome == "TRUE" and annotate == "TRUE":
    fin = open(vcf_in.replace(".vcf", "_genome_summary_formattedALL.txt"), "r")
elif exome == "TRUE" and annotate == "FALSE":
    fin = open(vcf_in, "r")

deNovo = open(final_out+"_deNovo.txt", "w")
rareHomozygous = open(final_out+"_homozygousRare.txt", "w")
compoundHet = open(final_out+"_compoundHet.txt", "w")
deNovoComplex = open(final_out+"_hemizygousGC.txt", "w")

counterDN = 0
counterRH = 0
counterCH = 0
counterHZGC = 0
counter = 0

#array and dicts for compound het finder
gene_found = []
gene_dict_pat = {}
gene_dict_mat = {}
wrote_dict = {}

header = fin.readline()
deNovo.write(header)
rareHomozygous.write(header)
compoundHet.write(header)
deNovoComplex.write(header)
headlist = header.replace("\n", "").split("\t")
backpopindex = vcfUtils.get_listindex(headlist, backpop)
offspring = headlist.index(args.offspringID)
father = headlist.index(args.fatherID)
mother = headlist.index(args.motherID)

MIEcounter = 0
n_confident = 0
oimputedcounter = 0
mimputedcounter = 0
fimputedcounter = 0

#read in and parse vcf formatted file
for line in fin.readlines():
    counter+=1
    linelist = line.split("\t")
    chrom = linelist[headlist.index("CHROM")]
    posit = linelist[headlist.index("POS")]
    info = linelist[headlist.index("INFO")]
    formatstats_avr_dp_rq = linelist[headlist.index("FORMAT")]
    #print "\npipeline info "+str(info)
    #print "\npipeline mq, mq0 "+str(mq)+" "+str(mq0)
    mq, mq0 = vcfUtils.parse_info(info)
    #print "\nstats field "+str(formatstats_avr_dp_rq)
    avr_num, dp_num, rq_num, gq_num = vcfUtils.parse_format(formatstats_avr_dp_rq)

#RQ	The phred scaled posterior probability that the sample is not identical to the reference.
#GQ	The standard VCF format genotype quality field. This is the phred scaled posterior score of the call. It is not necessarily the same as the QUAL column score.

    if ("#" in line) == 0:
        #print "\nline "+str(line)
        geno_father = vcfUtils.allele_coder(linelist[father], depth, altDepth, gtQual, pl, avr_num, dp_num, rq_num, gq_num, "RTG")
        geno_mother = vcfUtils.allele_coder(linelist[mother], depth, altDepth, gtQual, pl, avr_num, dp_num, rq_num, gq_num, "RTG")
        geno_offspring = vcfUtils.allele_coder(linelist[offspring], depth, altDepth, gtQual, pl, avr_num, dp_num, rq_num, gq_num, "RTG")
        alleles2 = geno_father+","+geno_mother+","+geno_offspring
        gene = linelist[headlist.index("Gene")]
        #This bit looks for imputed genotypes in the offspring only
        gto_list = linelist[offspring].split(":")
        if (int(gto_list[dp_num]) == 0):
            oimputedcounter+=1
        gto_list = linelist[mother].split(":")
        if (int(gto_list[dp_num]) == 0):
            mimputedcounter+=1
        gto_list = linelist[father].split(":")
        if (int(gto_list[dp_num]) == 0):
            fimputedcounter+=1
    if ((("3" in alleles2) == 0) and (("4" in alleles2) == 0)):
        tempcode = int(inmap[alleles2])
        n_confident+=1
        if tempcode in [1, 2, 5, 6, 8, 11, 15, 18, 20, 21, 24, 25]:
                MIEcounter+=1
    
    if ("#" in line) == 0 and vcfUtils.is_rare(line, freq, backpopindex) and mq >= mapQual and mq0 <= mapQual0:
        
        if ((("3" in alleles2) == 0) and (("4" in alleles2) == 0)):

            #simple de novo case
            if alleles2 == "0,0,1":
                deNovo.write(line)
                counterDN+=1

            #rare homozygous case
            elif alleles2 == "1,1,2":
                rareHomozygous.write(line)
                counterRH+=1

            #hemizygosity or gene conversion case
            elif alleles2 == "0,0,2" or alleles2 == "1,0,2" or alleles2 == "0,1,2" or alleles2 == "1,2,0" or alleles2 == "2,1,0":
                deNovoComplex.write(line)
                counterHZGC+=1

            #compound het case
            if gene_found.count(gene) >= 1:
                if gene_dict_mat.has_key(gene) and (alleles2 == "1,0,1"):
                    oldlist = gene_dict_mat[gene].split("\t")
		    chrom2 = oldlist[headlist.index("CHROM")]
                    posit2 = oldlist[headlist.index("POS")]
		    if wrote_dict.has_key(chrom2+":"+posit2) == 0:
                        compoundHet.write(gene_dict_mat[gene])
                        wrote_dict[chrom2+":"+posit2] = 1
                        counterCH+=1
                    compoundHet.write(line)
		    wrote_dict[chrom+":"+posit] = 1
                    gene_dict_pat[gene] = line
                    counterCH+=1
                elif gene_dict_pat.has_key(gene) and (alleles2 == "0,1,1"):
		    oldlist = gene_dict_pat[gene].split("\t")
		    chrom2 = oldlist[headlist.index("CHROM")]
		    posit2 = oldlist[headlist.index("POS")]
                    if wrote_dict.has_key(chrom2+":"+posit2) == 0:
                        compoundHet.write(gene_dict_pat[gene])
                        wrote_dict[chrom2+":"+posit2] = 1
                        counterCH+=1
                    compoundHet.write(line)
		    wrote_dict[chrom+":"+posit] = 1
                    gene_dict_mat[gene] = line
                    counterCH+=1
            elif (alleles2 == "0,1,1"):
                gene_dict_mat[gene] = line
                gene_found.append(gene)
            elif (alleles2 == "1,0,1"):
                gene_dict_pat[gene] = line
                gene_found.append(gene)

    if counter%500000 == 0:
        print "Processed "+str(counter)+" records"

print "\nRare de novo candidates: "+str(counterDN)
print "Rare homozygous candidates: "+str(counterRH)
print "Rare compound heterozygous candidates: "+str(counterCH)
print "Rare hemizygous or gene conversion candidates: "+str(counterHZGC)
print "MIE count: "+str(MIEcounter)
print "Confidently genotyped variants: "+str(n_confident)
print "Imputed variants in offspring: "+str(oimputedcounter)
print "Imputed variants in mother: "+str(mimputedcounter)
print "Imputed variants in father: "+str(fimputedcounter)

fin.close()
deNovo.close()
rareHomozygous.close()
compoundHet.close()
