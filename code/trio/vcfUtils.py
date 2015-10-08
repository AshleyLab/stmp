#Rick Dewey 5.22.12
#Last modified 8.22.13 James Priest
#module for general vcf utilities
#usage: utils module for vcf parsing
#!/usr/bin/env python

import os
import re
import sys

#gets proper head from vcf file, ignoring all other stuff
def get_head(infile):
    f1 = open(infile, "r")
    while True:
        line = f1.readline()
        if "CHROM" in line:
            head = line.replace("#", "")
            break
        if not line:
            print >> sys.stderr, "Error in vcfUtils.get_head - End of file reached prior to finding vcf header: vcf header is likely malformed"
            exit(1)
    f1.close()
    return head

#gets indexes from head of annotation file for downstream filtering/annotation
def get_listindex(head_list, viplist):
    indexlist = []
    tempindex = 0
    for item in head_list:
        if item in viplist:
            indexlist.append(tempindex)
        tempindex+=1
    return indexlist

#parses vcf format genotypes into allele codes, works for bialleleic positions, need to update for multiallelic positions, works for GATK standard vcf format only so far
def allele_coder(GT, d, alt, q, p, avr_n, dp_n, rq_n, gq_n, data_type):
    #deals with illumina format vcf files
    if (data_type == "ILLUMINA") or (data_type == "ILL"):
        likelihoodvec = gt_list[4].split(",")
        print "\nlike vec "+str(likelihoodvec)
        if ("./." in GT) == 0:
            gt_list = GT.split(":")
            gq = float(gt_list[3])
            temp = gt_list[1].split(",")
            if not ".:." in GT:
            #if int(temp[1]) != ".":
                al = int(temp[1])
                totalDepth = int(temp[0])+int(temp[1])
            # handles imputed genotypes from GATK phase by transmission tool which may impute genotypes without read information
            if ".:." in GT:
                al = 0
                totalDepth = 0
            
            if (gt_list[0] == "0/1" or gt_list[0] == "0|1" and totalDepth >= d and gq >= q and al >= alt and int(likelihoodvec[1]) <= p):
                alleles = "1"
            elif (gt_list[0] == "1|0" and totalDepth >= d and gq >= q and int(temp[0]) >= alt and int(likelihoodvec[0]) <= p):
                alleles = "1"
            elif (gt_list[0] == "1/1" or gt_list[0] == "1|1" and totalDepth >= d and gq >= q and al >= alt and int(likelihoodvec[2]) <= p):
                alleles = "2"
            elif (gt_list[0] == "0/0" or gt_list[0] == "0|0" and totalDepth >= d and gq >= q and int(likelihoodvec[0]) <= p):
                alleles = "0"

            #3 for alleles not meeting filtering criteria
            else:
                alleles = "3"

        #4 for alleles not called
        else:
            alleles = "4"
        return alleles

    #deals with complete genomics format vcf files
    if (data_type == "COMPLETE") or (data_type == "CG"):
        if ("./." in GT) == 0:
            gt_list = GT.split(":")
            if gt_list[0] == "0/1" or gt_list[0] == "0|1" or gt_list[0] == "1|0":
                alleles = "1"
            elif gt_list[0] == "1/1" or gt_list[0] == "1|1":
                alleles = "2"
            elif gt_list[0] == "0/0" or gt_list[0] == "0|0":
                alleles = "0"
            else:
                alleles = "3"

        #4 for alleles not called
        else:
            alleles = "4"
        return alleles
        
    #deals with RTG vcf files
    #AVR definition: adaptive variant rescoring value. It is a value between 0 and 1 that represents the probability that the variant for the sample is correct.
    #In the RTG version of the allele coder, likelihoodvec stores the AVR value, the --post_hoc_filter_cutoff represents the lower limit for inclusion
    #Cutoffs occur if the likelihoodvec is < p or pl
    if (data_type == "RTG") or (data_type == "rtg"):
        #print "\nGT "+str(GT)
        if ("./." in GT) == 0:
            gt_list = GT.split(":")
            likelihoodvec = float(gt_list[avr_n])
            gq = float(gt_list[gq_n])
            totalDepth = int(gt_list[dp_n])
            #print "\n p "+str(p)," AVR "+str(likelihoodvec)
            if (gt_list[0] == "1/1" or gt_list[0] == "1|1"):
                al = totalDepth * 0.5
            if totalDepth != 0:
                #In RTG formatted vcf files, the alternate allele depth not reported in the stats field
                al = int(99)
                #totalDepth = int(temp[0])+int(temp[1])
                #print "\ntotaldepth "+str(totalDepth), "\nal "+str(al), "\nlikelihoodvec "+str(likelihoodvec), "\ngq "+str(gq)
            # handles imputed genotypes from RTG which may impute genotypes without read information
            if totalDepth == 0:
                al = 0
                #print "\ntotaldepth "+str(totalDepth), "\nal "+str(al), "\nlikelihoodvec "+str(likelihoodvec), "\ngq "+str(gq)
            
            if (gt_list[0] == "0/1" or gt_list[0] == "0|1" and totalDepth >= d and gq >= q and al >= alt and likelihoodvec >= p):
                alleles = "1"
                #print "\ntotaldepth "+str(totalDepth), "\nal "+str(al), "\nlikelihoodvec "+str(likelihoodvec), "\ngq "+str(gq)
            elif (gt_list[0] == "1|0" and totalDepth >= d and gq >= q and al >= alt and likelihoodvec >= p):
                alleles = "1"
                #print "\ntotaldepth "+str(totalDepth), "\nal "+str(al), "\nlikelihoodvec "+ likelihoodvec), "\ngq "+str(gq)
            elif (gt_list[0] == "1/1" or gt_list[0] == "1|1" and totalDepth >= d and gq >= q and al >= alt and likelihoodvec >= p):
                alleles = "2"
                #print "\ntotaldepth "+str(totalDepth), "\nal "+str(al), "\nlikelihoodvec "+str(likelihoodvec), "\ngq "+str(gq)
            elif (gt_list[0] == "0/0" or gt_list[0] == "0|0" and totalDepth >= d and gq >= q and likelihoodvec >= p):
                #print "\ntotaldepth "+str(totalDepth), "\nal "+str(al), "\nlikelihoodvec "+str(likelihoodvec), "\ngq "+str(gq)
                alleles = "0"
                #print "\ntotaldepth "+str(totalDepth), "\nal "+str(al), "\nlikelihoodvec "+str(likelihoodvec), "\ngq "+str(gq)
            #3 for alleles not meeting filtering criteria
            else:
                alleles = "3"

        #4 for alleles not called
        else:
            alleles = "4"
        return alleles

    else:
        print >> sys.stderr, "Invalid input - unknown or unsupported data format: "+data_type
        exit(1)

	
#boolean, returns 1 for a variant with frequency below threshold (for all background populations consider) or novel or unknown allele frequency, 0 otherwise 
def is_rare(templine, freq_thresh, bp_indexlist):
    templinelist = templine.split("\t")
    rare_flag = 1
    for i in bp_indexlist:
        if templinelist[i] != "":
            if float(templinelist[i]) > float(freq_thresh):
                rare_flag = 0
    return rare_flag

#boolean, returns 1 for a variant with user-defined functional properties as annotated by annovar
def is_functional(templine, function_list):
    functional = 0
    if ("stoploss" in templine) and ("stoploss" in function_list):
        functional=1
    elif ("stopgain" in templine) and ("stopgain" in function_list):
        functional=1
    elif ("splicing" in templine) and ("splicing" in function_list):
        functional=1
    elif ("frameshift" in templine) and (("nonframeshift" in templine) == 0) and ("frameshift" in function_list):
        functional=1
    elif ("nonframeshift" in templine) and ("nonframeshift" in function_list):
        functional=1
    elif ("nonsynonymous" in templine) and ("nonsynonymous" in function_list):
        functional=1
    elif ("synonymous" in templine) and (("nonsynonymous" in templine) == 0) and ("synonymous" in function_list):
        functional=1
    elif ("exonic" in templine) and ("exonic" in function_list):
        functional=1
    elif ("intronic" in templine) and ("intronic" in function_list):
        functional=1
    elif ("UTR5" in templine) and ("UTR5" in function_list):
        functional=1
    elif ("UTR3" in templine) and ("UTR3" in function_list):
        functional=1
    elif ("ncRNA" in templine) and ("ncRNA" in function_list):
        functional=1
    elif ("upstream" in templine) and ("upstream" in function_list):
        functional=1
    elif ("intergenic" in templine) and ("intergenic" in function_list):
        functional=1
    elif ("downstream" in templine) and ("downstream" in function_list):
        functional=1
    return functional

#boolean, returns 1 for variant that is conserved according to user-defined criteria
def is_conserved(templine, cons_list, cons_criteria):
    templinelist = templine.split("\t")
    phyloP = templinelist[headlist.index("LJB_PhyloP_Pred")]
    gerp = templinelist[headlist.index("LJB_GERP++")]
    if gerp != "":
        gerp = float(gerp)
        if gerp >= 2.0 or phyloP == "C":
            return 1
        else:
            return 0
    elif phyloP == "C":
        return 1
    else:
        return 0

#boolean, returns 1 for variant that is pathogenic according to user-defined criteria
def is_pathogenic(templine, path_list, path_criteria):
    templinelist = templine.split("\t")
    pathogenic = 0
    sift = templinelist[headlist.index("LJB_SIFT_Pred")]
    pp2 = templinelist[headlist.index("LJB_PolyPhen2_Pred")]
    lrt = templinelist[headlist.index("LJB_LRT_Pred")]
    mt = templinelist[headlist.index("LJB_MutationTaster_Pred")]
    if sift == "D":
        pathogenic+=1
    if (pp2 == "P") or (pp2 == "D"):
        pathogenic+=1
    if lrt == "D":
        pathogenic+=1
    if (mt == "A") or (mt == "D"):
        pathogenic+=1
    if pathogenic >= int(nalg):
        return 1
    else:
        return 0
                
#parses info field of vcf file and returns tuple float for mq and mq0
def parse_info(infofield):
    infolist = infofield.split(";")
    #print "\nutils before infofield "+str(infofield)
    mapq = ""
    mapq0 = ""
    #print "\nutils before mapq, mapq0 "+str(mapq)+" "+str(mapq0)
    for element in infolist:
        if "MQ=" in element:
            mapq=float(element.replace("MQ=", ""))
        elif "MQ0=" in element:
            mapq0=float(element.replace("MQ0=", ""))
    #print "\nutils after infofield "+str(infofield)
    #print "\nutils mapq, mapq0 "+str(mapq)+" "+str(mapq0)
    return mapq, mapq0

#parses format field of vcf file to find location of AVR, DP, RQ
def parse_format(formatstats):
    temp2 = formatstats.split(':')
    counter = 0
    avr_num = 0 
    dp_num = 0
    rq_num = 0
    gq_num = 0
    for infos in temp2:
        if (str(infos) == "AVR"): 
            avr_num = int(counter)
        if (str(infos) == "DP"): 
            dp_num = int(counter)
        if (str(infos) == "RQ"): 
            rq_num = int(counter)
        if (str(infos) == "GQ"): 
            gq_num = int(counter)
        counter+=1
    return avr_num, dp_num, rq_num, gq_num
    


#merges vcf files split by chromosome, writing head from chromosome 1 to X only for now
def mergeFiles(fin_stem, f_out):
    chrom_arr = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]

    fout = open(f_out, "w")

    for chrom in chrom_arr:
        fin = open(fin_stem+"_chr"+chrom+"_annotated_hmmfiltered.txt", "r")
        head = fin.readline()
        if chrom == "1":
            fout.write(head)

        for line in fin.readlines():
            fout.write(line)
        fin.close()
        os.system("rm "+fin_stem+"_chr"+chrom+"_annotated_hmmfiltered.txt")
    fout.close()

#splits vcf file by chromosome, to X only for now
def splitFiles(f_in, fout_stem):
    chrom_arr = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]

    try:
        fin = open(f_in, "r")
    except IOError:
        print >> sys.stderr, "Error in vcfUtils.splitFiles: Cannot open input vcf file "+f_in
        exit(1)
                  
    head = get_head(f_in)
    try:
        f1 = open(fout_stem+"_chr1.txt", "w")
        f1.write(head)
        f2 = open(fout_stem+"_chr2.txt", "w")
        f2.write(head)
        f3 = open(fout_stem+"_chr3.txt", "w")
        f3.write(head)
        f4 = open(fout_stem+"_chr4.txt", "w")
        f4.write(head)
        f5 = open(fout_stem+"_chr5.txt", "w")
        f5.write(head)
        f6 = open(fout_stem+"_chr6.txt", "w")
        f6.write(head)
        f7 = open(fout_stem+"_chr7.txt", "w")
        f7.write(head)
        f8 = open(fout_stem+"_chr8.txt", "w")
        f8.write(head)
        f9 = open(fout_stem+"_chr9.txt", "w")
        f9.write(head)
        f10 = open(fout_stem+"_chr10.txt", "w")
        f10.write(head)
        f11 = open(fout_stem+"_chr11.txt", "w")
        f11.write(head)
        f12 = open(fout_stem+"_chr12.txt", "w")
        f12.write(head)
        f13 = open(fout_stem+"_chr13.txt", "w")
        f13.write(head)
        f14 = open(fout_stem+"_chr14.txt", "w")
        f14.write(head)
        f15 = open(fout_stem+"_chr15.txt", "w")
        f15.write(head)
        f16 = open(fout_stem+"_chr16.txt", "w")
        f16.write(head)
        f17 = open(fout_stem+"_chr17.txt", "w")
        f17.write(head)
        f18 = open(fout_stem+"_chr18.txt", "w")
        f18.write(head)
        f19 = open(fout_stem+"_chr19.txt", "w")
        f19.write(head)
        f20 = open(fout_stem+"_chr20.txt", "w")
        f20.write(head)
        f21 = open(fout_stem+"_chr21.txt", "w")
        f21.write(head)
        f22 = open(fout_stem+"_chr22.txt", "w")
        f22.write(head)
        fX = open(fout_stem+"_chrX.txt", "w")
        fX.write(head)
    except IOError:
        print >> sys.stderr, "Error in vcfUtils.splitFiles: Improper output file specification "+f_in 
        exit(1)
                  
    while 1:
        temp = fin.readline()
        if not temp:
            break
        else:
            if ("chr1\t" in temp) and (('#' in temp) == 0):
                f1.write(temp)
            if ("chr2\t" in temp) and (('#' in temp) == 0):
                f2.write(temp)
            if ("chr3\t" in temp) and (('#' in temp) == 0):
                f3.write(temp)
            if ("chr4\t" in temp) and (('#' in temp) == 0):
                f4.write(temp)
            if ("chr5\t" in temp) and (('#' in temp) == 0):
                f5.write(temp)
            if ("chr6\t" in temp) and (('#' in temp) == 0):
                f6.write(temp)
            if ("chr7\t" in temp) and (('#' in temp) == 0):
                f7.write(temp)
            if ("chr8\t" in temp) and (('#' in temp) == 0):
                f8.write(temp)
            if ("chr9\t" in temp) and (('#' in temp) == 0):
                f9.write(temp)
            if ("chr10\t" in temp) and (('#' in temp) == 0):
                f10.write(temp)
            if ("chr11\t" in temp) and (('#' in temp) == 0):
                f11.write(temp)
            if ("chr12\t" in temp) and (('#' in temp) == 0):
                f12.write(temp)
            if ("chr13\t" in temp) and (('#' in temp) == 0):
                f13.write(temp)
            if ("chr14\t" in temp) and (('#' in temp) == 0):
                f14.write(temp)
            if ("chr15\t" in temp) and (('#' in temp) == 0):
                f15.write(temp)
            if ("chr16\t" in temp) and (('#' in temp) == 0):
                f16.write(temp)
            if ("chr17\t" in temp) and (('#' in temp) == 0):
                f17.write(temp)
            if ("chr18\t" in temp) and (('#' in temp) == 0):
                f18.write(temp)
            if ("chr19\t" in temp) and (('#' in temp) == 0):
                f19.write(temp)
            if ("chr20\t" in temp) and (('#' in temp) == 0):
                f20.write(temp)
            if ("chr21\t" in temp) and (('#' in temp) == 0):
                f21.write(temp)
            if ("chr22\t" in temp) and (('#' in temp) == 0):
                f22.write(temp)
            if ("chrX\t" in temp) and (('#' in temp) == 0):
                fX.write(temp)
