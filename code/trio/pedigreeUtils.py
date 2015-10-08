#Rick Dewey 6.28.12
#Last modified 8.26.13 James Priest
#module for general pedigree utilities
#usage: utils module for general pedigree utilities
#!/usr/bin/env python

import os
import re
import sys
import math
import vcfUtils
#import hmmUtils
import time
import numpy as np


#traverses .ped files and returns founders, a dictionary with parents and number of offspring (trios), all trios, quartets, and parents with one child only (unique_trios), dictionary with affected status
def traverse(pfile):
    fam_arr = []
    parent_arr = []
    parent_dict = {}
    last_trio = {}
    affected = {}
    
    trios = []
    unique_trios = []
    quartets = []
    founders = []

    try:
        pedfile = open(pfile, "r")
    except IOError:
        print >> sys.stderr, "Error in pedigreeUtils.traverse: Invalid pedigree file"
        exit(1)
        
    header = pedfile.readline()
    
    for line in pedfile:
        linelist = line.replace("\n", "").split("\t")
        family_id = linelist[0]
        child = linelist[1]
        father = linelist[2]
        mother = linelist[3]
        gender = linelist[4]
        affected[child] = linelist[5]

        if (father == "--") and (mother == "--"):
            founders.append(child)
        else:
            
            #trios sharing same parent are identified as quartets
            if ((father+":"+mother) in parent_arr):
                parent_dict[father+":"+mother]+=1
                trios.append(father+":"+mother+":"+child)
                quartets.append(father+":"+mother+":"+last_trio[father+":"+mother]+":"+child)
                last_trio[father+":"+mother] = child

                #remove trios that form quartets from unique_trios
                for element in unique_trios:
                    if father+":"+mother in element:
                        unique_trios.remove(element)
            else:
                parent_arr.append(father+":"+mother)
                parent_dict[father+":"+mother] = 1
                trios.append(father+":"+mother+":"+child)
                last_trio[father+":"+mother] = child
                unique_trios.append(father+":"+mother+":"+child)
    
    return founders, parent_dict, trios, quartets, unique_trios, affected

#creates inheritance states for quartets: 0 = haploidentical paternal, 1 = haploidentical maternal, 2 = identical, 3 = nonidentical; and two error states: 4 = compression, 5 = mie rich
def phase_quartet(fatherID, motherID, offspring1ID, offspring2ID, vcf_in_stem, depth, altDepth, gtQual, mapQual, mapQual0, pl, mat_path, error, compression, mie, selfprob, platform):
    nonselfprob = str((1-float(selfprob))/float(5))

    #create matrices for HMM, all probabilities in log10 space to avoid underflow from floating point
    transmat, emitmat, startp, states, inmap = hmmUtils.make_hmm(mat_path, selfprob, compression, mie, error, "q")

    #temp file for inheritance state
    phase_out = open(".".join([vcf_in_stem, fatherID, motherID, offspring1ID, offspring2ID, "phase_out.txt"]), "w")

    #chrom_arr = ["21", "22"]
    chrom_arr = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]

    start_time=time.clock()
    for chromh in chrom_arr:
        print "\nWorking on chromosome "+str(chromh)+"\n"
        fin = open(vcf_in_stem+"_chr"+chromh+".txt", "r")
    
        header = fin.readline()

        #grab proper columns from header for variants
        headlist = header.replace("\n", "").split("\t")
        try:
            offspring1 = headlist.index(offspring1ID)
            offspring2 = headlist.index(offspring2ID)
            father = headlist.index(fatherID)
            mother = headlist.index(motherID)
        except ValueError:
            return "Error in pedigreeUtils.phase_quartet: Invalid input for quartet "+":".join([fatherID, motherID, offspring1ID, offspring2ID])+" - Specified subject ID does not match vcf file header, consider malformed header or pedigree file"
            
        #arrays for observations and positions
        obs_seq = []
        position_clean = []
        position_all = []
        res_dict = {}

        #counters for error modes (DN = de novo allele, HZGC = apparent gene conversions/hemizygosity
        counterDN1 = 0
        counterHZGC1 = 0
        counterDN2 = 0
        counterHZGC2 = 0
        counter = 0

        #counters for number of informative allele assortments for determing paternity/maternity
        inform_parents = 0
        non_paternal1 = 0
        non_maternal1 = 0
        non_paternal2 = 0
        non_maternal2 = 0

        for line in fin.readlines():
            if ("#" in line) == 0:
                linelist = line.replace("\n", "").split("\t")
                chrom = linelist[headlist.index("CHROM")]
                posit = linelist[headlist.index("POS")]
                info = linelist[headlist.index("INFO")]
                position_all.append(posit)

                #temp variable to flag if line should be used, and also check for unsupported data format (currently supports CG and Illumina GATK-style)
                goline = 0
                if ((platform == "ILLUMINA") or (platform == "ILL")):
                    try:
                        mq, mq0 = vcfUtils.parse_info(info)
                        if mq >= mapQual and mq0 <= mapQual0 and ("PASS" in line):
                            goline = 1
                    except UnboundLocalError:
                        return "Error in pedigreeUtils.phase_quartet for quartet "+":".join([fatherID, motherID, offspring1ID, offspring2ID])+" - Invalid input format: Data does not appear to be in GATK style from Illumina platform"
                elif (platform == "COMPLETE") or (platform == "CG"):
                    goline = 1
                elif (platform == "RTG") or (platform == "rtg"):
                    goline = 1
                elif (platform != "COMPLETE") and (platform != "CG") and (platform != "ILLUMINA") and (platform != "ILL"):
                    return "Error in pedigreeUtils.phase_quartet for quartet "+":".join([fatherID, motherID, offspring1ID, offspring2ID])+": Invalid input - unknown or unsupported data format: "+platform

                if goline == 1:
                    geno_father = vcfUtils.allele_coder(linelist[father], depth, altDepth, gtQual, pl, platform)
                    geno_mother = vcfUtils.allele_coder(linelist[mother], depth, altDepth, gtQual, pl, platform)
                    geno_offspring1 = vcfUtils.allele_coder(linelist[offspring1], depth, altDepth, gtQual, pl, platform)
                    geno_offspring2 = vcfUtils.allele_coder(linelist[offspring2], depth, altDepth, gtQual, pl, platform)

                    alleles2 = geno_father+","+geno_mother+","+geno_offspring1+","+geno_offspring2
                    if ((("3" in alleles2) == 0) and (("4" in alleles2) == 0)):
                        counter+=1
                        
                        #checks for paternity, maternity
                        if geno_father == "0" and geno_mother == "2":
                            inform_parents += 1
                            if geno_offspring1 == "2":
                                non_paternal1 += 1
                            if geno_offspring1 == "0":
                                non_maternal1 += 1
                            if geno_offspring2 == "2":
                                non_paternal2 += 1
                            if geno_offspring2 == "0":
                                non_maternal2 += 1

                        elif geno_father == "2" and geno_mother == "0":
                            inform_parents += 1
                            if geno_offspring1 == "0":
                                non_paternal1 += 1
                            if geno_offspring1 == "2":
                                non_maternal1 += 1
                            if geno_offspring2 == "0":
                                non_paternal2 += 1
                            if geno_offspring2 == "2":
                                non_maternal2 += 1
                    
                        #simple de novo case
                        if alleles2 == "0,0,1,1":
                            counterDN1+=1
                            counterDN2+=1
                        elif alleles2 == "0,0,1,0":
                            counterDN1+=1
                        elif alleles2 == "0,0,0,1":
                            counterDN2+=1
                        
                        #hemizygosity or gene conversion case
                        elif (geno_father == "0" and geno_mother == "0" and geno_offspring1 == "2") or (geno_father == "1" and geno_mother == "0" and geno_offspring1 == "2") or (geno_father == "0" and geno_mother == "1" and geno_offspring1 == "2") or (geno_father == "1" and geno_mother == "2" and geno_offspring1 == "0") or (geno_father == "2" and geno_mother == "1" and geno_offspring1 == "0"):
                            counterHZGC1+=1
                        elif (geno_father == "0" and geno_mother == "0" and geno_offspring2 == "2") or (geno_father == "1" and geno_mother == "0" and geno_offspring2 == "2") or (geno_father == "0" and geno_mother == "1" and geno_offspring2 == "2") or (geno_father == "1" and geno_mother == "2" and geno_offspring2 == "0") or (geno_father == "2" and geno_mother == "1" and geno_offspring2 == "0"):
                            counterHZGC2+=1

                        #creates an array of clean observations for HMM
                        if inmap.has_key(alleles2):
                            obs_seq.append(int(inmap[alleles2]))
                            position_clean.append(posit)

        fin.close()
                
        cent_np1 = float(non_paternal1)/float(inform_parents)
        cent_nm1 = float(non_maternal1)/float(inform_parents)
        cent_DN1 = float(counterDN1+counterHZGC1)/float(counter)

        cent_np2 = float(non_paternal2)/float(inform_parents)
        cent_nm2 = float(non_maternal2)/float(inform_parents)
        cent_DN2 = float(counterDN2+counterHZGC2)/float(counter)

        print "\nAllele QC summary for quartet: "+str(":".join([fatherID, motherID, offspring1ID, offspring2ID])+"\n")
        print "\nNon-transmitted paternal alleles for child 1 (%): "+str(cent_np1*float(100))
        print "Non-transmitted maternal alleles for child 1 (%): "+str(cent_nm1*float(100))
        print "Novel alleles for child 1 (%): "+str(cent_DN1*float(100))
        print "\nNon-transmitted paternal alleles for child 2 (%): "+str(cent_np2*float(100))
        print "Non-transmitted maternal alleles for child 2 (%): "+str(cent_nm2*float(100))
        print "Novel alleles for child 2 (%): "+str(cent_DN2*float(100))

        #warning for high rates of non-transmission, does not kill process
        if (cent_np1) > 0.05 and (chromh != 'X'):
            print "\nWARNING! High rate of non-transmission of paternal alleles in child "+offspring1ID+": consider sample mix-up or non-paternity...\n"
        if (cent_np2) > 0.05 and (chromh != 'X'):
            print "\nWARNING! High rate of non-transmission of paternal alleles in child "+offspring2ID+": consider sample mix-up or non-paternity...\n"
        if (cent_nm1) > 0.05 and (chromh != 'X'):
            print "\nWARNING! High rate of non-transmission of maternal alleles in child "+offspring1ID+": consider sample mix-up or non-maternity...\n"
        if (cent_nm2) > 0.05 and (chromh != 'X'):
            print "\nWARNING! High rate of non-transmission of maternal alleles in child "+offspring2ID+": consider sample mix-up or non-maternity...\n"
        if (cent_DN1) > 0.05 and (chromh != 'X'):
            print "\nWARNING! High rate of novel alleles in child "+offspring1ID+": consider sample mix-up or non-maternity or non-paternity...\n"
        if (cent_DN2) > 0.05 and (chromh != 'X'):
            print "\nWARNING! High rate of novel alleles in child "+offspring2ID+": consider sample mix-up or non-maternity or non-paternity...\n"

        #run Viterbi algorithm on HMM defined above
        vit = hmmUtils.viterbi(np.array(obs_seq), states, startp, transmat, emitmat)
        vit_path = vit[1]

        #bind viterbi path to position dict for annotation
        for i in range(0, len(vit_path)):
            res_dict[position_clean[i]] = vit_path[i]

        #temp variables to allow bridging alleles not used in HMM
        firstpath = lastpath = vit_path[0]

        #create output dict with state for every position, including those not used in HMM
        for i in position_all:
            if res_dict.has_key(i):
                phase_out.write("chr"+chromh+":"+str(i)+"\t"+str(res_dict[i])+"\n")
                lastpath = res_dict[i]
            else:
                phase_out.write("chr"+chromh+":"+str(i)+"\t"+str(lastpath)+"\n")
    end_time = time.clock()
    phase_out.close()
    return "Normal exit status for quartet "+":".join([fatherID, motherID, offspring1ID, offspring2ID])+". Processing time: "+str((float(end_time)-float(start_time))/float(60))+" minutes"

#wrapper for QA hmm for trios: 0 = mierich, 1 = compression, 2 = good data
def phase_trio(fatherID, motherID, offspringID, vcf_in_stem, depth, altDepth, gtQual, mapQual, mapQual0, pl, mat_path, error, compression, mie, selfprob, platform):

    nonselfprob = str((float(1)-float(selfprob))/float(2))

    #create matrices for HMM, all probabilities in log10 space to avoid underflow
    transmat, emitmat, startp, states, inmap = hmmUtils.make_hmm(mat_path, selfprob, compression, mie, error, "t")

    #chrom_arr = ["21", "22"]
    chrom_arr = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]

    #temp file for inheritance state
    phase_out = open(".".join([vcf_in_stem, fatherID, motherID, offspringID, "phase_out.txt"]), "w")
    start_time=time.clock()
    
    for chromh in chrom_arr:
        fin = open(vcf_in_stem+"_chr"+chromh+".txt", "r")
        print "\nWorking on chromosome "+str(chromh)+"\n"
        header = fin.readline()

        #grab proper columns from header for variants, need error message for improper input
        headlist = header.replace("\n", "").split("\t")
        try: 
            offspring = headlist.index(offspringID)
            father = headlist.index(fatherID)
            mother = headlist.index(motherID)
        except ValueError:
            return "Error in pedigreeUtils.phase_trio for trio "+":".join([fatherID, motherID, offspringID])+": Invalid input - Specified subject ID does not match vcf file header, consider malformed header or pedigree file"
    
        obs_seq = []
        position_clean = []
        position_all = []
        res_dict = {}

        counterDN = 0
        counterHZGC = 0
        counter = 0

        #counter for number of informative allele assortments for determing paternity/maternity
        inform_parents = 0
        non_paternal = 0
        non_maternal = 0

        for line in fin.readlines():
            if ("#" in line) == 0:
                counter+=1
                linelist = line.replace("\n", "").split("\t")
                chrom = linelist[headlist.index("CHROM")]
                posit = linelist[headlist.index("POS")]
                info = linelist[headlist.index("INFO")]
                position_all.append(posit)
    
                #temp variable to flag if line should be used, and also check for unsupported data format (currently supports CG and Illumina GATK-style)
                goline = 0
                if ((platform == "ILLUMINA") or (platform == "ILL")):
                    try:
                        mq, mq0 = vcfUtils.parse_info(info)
                        if mq >= mapQual and mq0 <= mapQual0 and ("PASS" in line):
                            goline = 1
                    except UnboundLocalError:
                        return "Error in pedigreeUtils.phase_trio for trio "+":".join([fatherID, motherID, offspringID])+": Data does not appear to be in GATK style from Illumina platform"
                elif (platform == "COMPLETE") or (platform == "CG"):
                    goline = 1
                elif (platform == "RTG") or (platform == "rtg"):
                    goline = 1
                elif (platform != "COMPLETE") and (platform != "CG") and (platform != "ILLUMINA") and (platform != "ILL"):
                    return "Error in pedigreeUtils.phase_trio for trio "+":".join([fatherID, motherID, offspringID])+": Invalid input - unknown or unsupported data format: "+platform

                if goline == 1:
        
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
                            position_clean.append(posit)

        fin.close()
                
        cent_np = float(non_paternal)/float(inform_parents)
        cent_nm = float(non_maternal)/float(inform_parents)
        cent_DN = float(counterDN+counterHZGC)/float(counter)

        print "\nAllele QC summary for trio: "+str(":".join([fatherID, motherID, offspringID])+"\n")
        print "\nNon-transmitted paternal alleles (%): "+str(cent_np*float(100))
        print "Non-transmitted maternal alleles (%): "+str(cent_nm*float(100))
        print "Novel alleles (%): "+str(cent_DN*float(100))

        #warning for high rates of non-transmission on autosomes, does not kill process
        if (cent_np > 0.05) and (chromh != 'X'):
            print "\nWARNING! High rate of non-transmission of paternal alleles: consider sample mix-up or non-paternity...\n"
        if (cent_nm > 0.05) and (chromh != 'X'):
            print "\nWARNING! High rate of non-transmission of maternal alleles: consider sample mix-up or non-maternity...\n"
        if (cent_DN > 0.05) and (chromh != 'X'):
            print "\nWARNING! High rate of novel alleles in child: consider sample mix-up or non-maternity or non-paternity...\n"

        fin.close()

        #run Viterbi algorithm on HMM defined above
        vit = hmmUtils.viterbi(np.array(obs_seq), states, startp, transmat, emitmat)
        vit_path = vit[1]

        #bind viterbi path to position dict for annotation
        for i in range(0, len(vit_path)):
            res_dict[position_clean[i]] = vit_path[i]

        #temp variables to allow bridging alleles not used in HMM
        firstpath = lastpath = vit_path[0]

        #create output dict with state for every position, including those not used in HMM
        for i in position_all:
            if res_dict.has_key(i):
                phase_out.write("chr"+chromh+":"+str(i)+"\t"+str(res_dict[i])+"\n")
                lastpath = res_dict[i]
            else:
                phase_out.write("chr"+chromh+":"+str(i)+"\t"+str(lastpath)+"\n")
    end_time = time.clock()
    phase_out.close()
    return "Normal exit status for trio "+":".join([fatherID, motherID, offspringID])+". Processing time: "+str((float(end_time)-float(start_time))/float(60))+" minutes"
