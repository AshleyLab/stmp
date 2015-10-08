#Rick Dewey 5.22.12
#Last modified 1.26.14 Rick Dewey
#module for general vcf utilities
#usage: utils module for vcf parsing
#!/usr/bin/env python

import os
import re
import sys
import vcfHeaders

#gets proper head from vcf file, ignoring all other stuff
def get_head(infile):
    f1 = open(infile, "r")
    while True:
        line = f1.readline()
        if "#CHROM" in line:
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

#parses vcf format genotypes into allele codes, works for bialleleic positions, need to update for multiallelic positions
def allele_coder(GT, d, alt, q, p, data_type):

    #deals with illumina format vcf files as processed by gatk
    if (data_type == "ILLUMINA") or (data_type == "ILL"):
        if ("./." in GT) == 0:
            gt_list = GT.split(":")
            gq = float(gt_list[3])
            temp = gt_list[1].split(",")
            al = int(temp[1])
            totalDepth = int(temp[0])+int(temp[1])
            likelihoodvec = gt_list[4].split(",")

            if (gt_list[0] == "0/1" and totalDepth >= d and gq >= q and al >= alt and int(likelihoodvec[1]) <= p):
                alleles = "1"
            elif (gt_list[0] == "1/1" and totalDepth >= d and gq >= q and al >= alt and int(likelihoodvec[2]) <= p):
                alleles = "2"
	    elif (gt_list[0] == "0/0" and totalDepth >= d and gq >= q and int(likelihoodvec[0]) <= p):
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
            if gt_list[0] == "0/1":
                alleles = "1"
            elif gt_list[0] == "1/1":
                alleles = "2"
	    elif gt_list[0] == "0/0":
                alleles = "0"
            elif gt_list[0] == "0/.":
                alleles = "0"
            elif gt_list[0] == "1/.":
                alleles = "2"
            else:
                alleles = "3"

        #4 for alleles not called
        else:
            alleles = "4"
        return alleles
    else:
        print >> sys.stderr, "Invalid input - unknown or unsupported data format: "+data_type
        exit(1)

    #deals with real time genomics format vcf files

#boolean, returns 1 for a variant with frequency below threshold (for all background populations consider) or novel or unknown allele frequency, 0 otherwise 
def is_rare(templine, freq_thresh, bp_indexlist):
    templinelist = templine.split("\t")
    rare_flag = 1 # rare
    for i in bp_indexlist:
        if templinelist[i] != "":
            # should throw an exception if templinelist[i] isn't a float
            if float(templinelist[i]) > float(freq_thresh):
                rare_flag = 0 # not rare
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
def is_conserved(templine, headlist, thresh):
    total = 0
    templinelist = templine.split("\t")
    phyloP = ''
    if(vcfHeaders.kPhyloP_pred in headlist):
        phyloP = templinelist[headlist.index(vcfHeaders.kPhyloP_pred)]
    else:
        print 'warning: ' + vcfHeaders.kPhyloP_pred + ' not found in annotated header'
        
    gerp = ''
    if(vcfHeaders.kGerp in headlist):
        gerp = templinelist[headlist.index(vcfHeaders.kGerp)]
    else:
        print 'warning: ' + vcfHeaders.kGerp + ' not found in annotated header'
        
    phast_cons = ''
    if(vcfHeaders.kPhastConsElts46Way in headlist):
        phast_cons = templinelist[headlist.index(vcfHeaders.kPhastConsElts46Way)]
    else:
        print 'warning: ' + vcfHeaders.kPhastConsElts46Way + ' not found in annotated header'
    
    wg_gerp = ''
    if(vcfHeaders.kWg_gerp in headlist): 
        wg_gerp = templinelist[headlist.index(vcfHeaders.kWg_gerp)]
    else:
        print 'warning: ' + vcfHeaders.kWg_gerp + ' not found in annotated header'

    if gerp != "" and phyloP != "":
        gerp = float(gerp)
        if gerp >= 2.0:
            total+=1
    if phyloP == "C":
        total+=1
    if phast_cons!= "":
        if int(phast_cons) >= 250:
            total+=1
    if wg_gerp != "":
        if float(wg_gerp) >= 2.0:
            total+=1
    if total >= thresh:
        return 1
    else:
        return 0


def getClinvarInfoCol(templine, headlist):
    templinelist = templine.split("\t")
    clinvarInfo = templinelist[headlist.index(vcfHeaders.kClinvar)] if vcfHeaders.kClinvar in headlist else ''
    if(vcfHeaders.kClinvar not in headlist):
        print 'warning: ' + vcfHeaders.kClinvar + ' not found in annotated header'
    return clinvarInfo

# TODO be consistent -- when doing clinvar star annotation, need clinvar clinrevstatus without "clinvar_" prefix, but when doing tiering with final annotated VCF, need clinvar clinsig with "clinvar_" prefix
def getClinvarClinsigHeader(yaml_cmds):
    return yaml_cmds['clinvar']['Annotation']+'_'+vcfHeaders.kClinvarClinSigHeader
#     return vcfHeaders.kClinvarClinSigHeader
def getClinvarClinRevStatusHeader(yaml_cmds):
#     return yaml_cmds['clinvar']['Annotation']+'_'+vcfHeaders.kClinvarClinRevStatusHeader
    return vcfHeaders.kClinvarClinRevStatusHeader

# return clinsig value(s) as a string (may have multiple reported values in clinvar)
def getClinvarClinsigStr(templine, headlist, yaml_cmds):
    clinvarClinsigHeader = getClinvarClinsigHeader(yaml_cmds)
    if(clinvarClinsigHeader in headlist):
        idx = headlist.index(clinvarClinsigHeader)
        tempLineElts = templine.rstrip("\n").split("\t")
        return tempLineElts[idx]
    else:
        raise IndexError('could not get index of clinvar clinical significance column (' + clinvarClinsigHeader + ') in annotated file')
        print 'error could not get index of clinvar clinical review status column in annotated file'
        return ''
#     clinvarInfo = getClinvarInfoCol(templine, headlist)
#     if('CLNSIG' in clinvarInfo):
#         None # TODO finish
#     # else
#     return ''
def clinvarClinsigStr2Array(clinsig_str):
    return clinsig_str.split('|')
def getClinvarClinsigArray(templine, headlist):
    return clinvarClinsigStr2Array(getClinvarClinsigStr(templine, headlist))

def isClinvarPathogenicOrLikelyPathogenic(line, headlist, yaml_cmds):
    clinsigStr = getClinvarClinsigStr(line, headlist, yaml_cmds)
    if(str(vcfHeaders.kCLINVAR_CLINSIG_LIKELY_PATHOGENIC_CODE) in clinsigStr or str(vcfHeaders.kCLINVAR_CLINSIG_PATHOGENIC_CODE) in clinsigStr):
        return True
    return False

#return clinrevstatus values as a string (similar to clinsig above)
def getClinvarClinReviewStatusStr(templine, headlist, yaml_cmds):
    clinvarClinRevStatusStr = getClinvarClinRevStatusHeader(yaml_cmds)
    if(clinvarClinRevStatusStr in headlist):
        idx = headlist.index(clinvarClinRevStatusStr)
        tempLineElts = templine.rstrip("\n").split("\t")
        return tempLineElts[idx]
    else:
        raise IndexError('could not get index of clinvar clinical review status column (' + clinvarClinRevStatusStr +') in annotated file')
        print 'error could not get index of clinvar clinical review status column in annotated file'
        return ''
# clinrevstatus convenience methods
def clinvarClinRevStatusStr2Array(clinrevstat_str):
    return clinrevstat_str.split('|')
def getClinvarClinReviewStatusArray(templine, headlist):
    return clinvarClinRevStatusStr2Array(getClinvarClinReviewStatusStr(templine, headlist))


# computes clinvar stars based on clinsig + clin review status (based on revised ClinVar mid 2015 criteria)
# from http://www.ncbi.nlm.nih.gov/clinvar/docs/variation_report/#review_status
# 0 stars: No submitter provided an interpretation with assertion criteria (no assertion criteria provided), or no interpretation was provided (no assertion provided)
# 1 star: At least one submitter provided an interpretation with assertion criteria (criteria provided, single submitter) or multiple submitters provided assertion criteria but there are conflicting interpretations in which case the independent values are enumerated for clinical significance (criteria provided, conflicting interpretations)
# 2 stars: Two or more submitters providing assertion criteria provided the same interpretation (criteria provided, multiple submitters, no conflicts)
# 3 stars: reviewed by expert panel (http://www.ncbi.nlm.nih.gov/clinvar/docs/review_guidelines/)
# 4 stars: practice guideline (http://www.ncbi.nlm.nih.gov/clinvar/docs/review_guidelines/)
def clinvarStars(templine, headlist, yaml_cmds):
    clinRevStatStr = getClinvarClinReviewStatusStr(templine, headlist, yaml_cmds)
    
    # 4 star: practice guideline
    if(clinRevStatStr.find('practice_guideline') != -1):
        return 4
    
    # 3 star: expert panel
    if(clinRevStatStr.find('exp') != -1):
        return 3
    
    #0, 1, and 2 star: check assertion criteria and whether assertions agree
    numAssertionCriteriaProvided = clinRevStatStr.count('criteria_provided') - clinRevStatStr.count('no_assertion_criteria_provided')
    if(numAssertionCriteriaProvided == 0):
        return 0
    if(numAssertionCriteriaProvided == 1 and clinRevStatStr.find('criteria_provided\x2c_multiple_submitters') == -1):
        return 1
    if(numAssertionCriteriaProvided > 1 or clinRevStatStr.find('criteria_provided\x2c_multiple_submitters') != -1):
        if clinRevStatStr.find('conflicting_interpretations') != -1:
            return 1
        #else, no conflicting interpretations
        return 2


#boolean, returns 1 for variant that is pathogenic according to user-defined criteria
def is_pathogenic(templine, headlist, nalg):
    templinelist = templine.split("\t")
    pathogenic = 0
    
    sift = templinelist[headlist.index(vcfHeaders.kSiftPred)] if vcfHeaders.kSiftPred in headlist else ''
    pp2 = templinelist[headlist.index(vcfHeaders.kPolyphen2Pred)] if vcfHeaders.kPolyphen2Pred in headlist else ''
    lrt = templinelist[headlist.index(vcfHeaders.kLrtPred)] if vcfHeaders.kLrtPred in headlist else ''
    mt = templinelist[headlist.index(vcfHeaders.kMutationTasterPred)] if vcfHeaders.kMutationTasterPred in headlist else ''
    pp2_2 = templinelist[headlist.index(vcfHeaders.kPolyphen2Pred_2)] if vcfHeaders.kPolyphen2Pred_2 in headlist else ''
    
    # check for header mismatches
    if(vcfHeaders.kSiftPred not in headlist):
        print 'warning: ' + vcfHeaders.kSiftPred + ' not found in annotated header'
    if(vcfHeaders.kPolyphen2Pred not in headlist):
        print 'warning: ' + vcfHeaders.kPolyphen2Pred + ' not found in annotated header'
    if(vcfHeaders.kLrtPred not in headlist):
        print 'warning: ' + vcfHeaders.kLrtPred + ' not found in annotated header'
    if(vcfHeaders.kMutationTasterPred not in headlist):
        print 'warning: ' + vcfHeaders.kMutationTasterPred + ' not found in annotated header'
    if(vcfHeaders.kPolyphen2Pred_2 not in headlist):
        print 'warning: ' + vcfHeaders.kPolyphen2Pred_2 + ' not found in annotated header (though currently unused)'
    
    if sift == "D":
        pathogenic+=1
        #debug
        print 'SIFT=D'
    if (pp2 == "P") or (pp2 == "D"):
        pathogenic+=1
        #debug 
        print 'pp=P or D'
    if lrt == "D":
        pathogenic+=1
        # debug
        print 'lrt=D'
    if (mt == "A") or (mt == "D"):
        pathogenic+=1
        #debug
        print 'mt=A or D'
    if pathogenic >= int(nalg):
        #debug
        print 'is pathogenic'
        return 1
    else:
        print 'not pathogenic'
        return 0


#CURRENTLY UNUSED
#finds variants with regulatory annotations
def is_regulatory(templine, headlist):
    templinelist = templine.split("\t")
    regulatory = 0
    dnase = 0
    adult_enhancer = 0
    fetal_enhancer = 0
    tfbs = 0
    transfac=0
    mirna_coding=0
    mirna_target=0
    if templinelist[headlist.index(vcfHeaders.kEncode_dnaseCluster)] != "":
        if int(templinelist[headlist.index(vcfHeaders.kEncode_dnaseCluster)]) >= 300:
            regulatory=1
            dnase = 1
    if templinelist[headlist.index("encode_tfbsChip_score")] != "": # warning not in annotation
        if int(templinelist[headlist.index("encode_tfbsChip_score")]) >= 300:
            regulatory=1
            tfbs=1
    if templinelist[headlist.index("heart_adult_enhancer")] != "": # warning not in annotation
        regulatory=1
        adult_enhancer=1
    if templinelist[headlist.index("heart_fetal_enhancer")] != "": # warning not in annotation
        regulatory=1
        fetal_enhancer=1
    if templinelist[headlist.index("transFac_score")] != "": # warning not in annotation
        regulatory=1
        transfac=1
    if templinelist[headlist.index("miRNA")] != "": # warning not in annotation
        regulatory=1
        mirna_coding=1
    if templinelist[headlist.index("targetScan")] != "": # warning not in annotation
        regulatory=1
        mirna_target=1
    return regulatory, dnase, tfbs, adult_enhancer, fetal_enhancer, transfac, mirna_coding, mirna_target
    

# CURRENTLY UNUSED
#finds variants in topologically central location in network according to user defined criteria
def is_central_mod(templine, headlist, rank_thresh, phenotype):
    templinelist = templine.split("\t")
    if phenotype == "normal":
        if templinelist[headlist.index("normal_heart_Kin")] != "": # warning not in annotation
            if float(templinelist[headlist.index("normal_heart_Kin")]) >= rank_thresh:
                return 1
            else:
                return 0
        else:
            return 0
    if phenotype == "hypertrophy":
        if templinelist[headlist.index("hypertrophy_heart_Kin")] != "": # warning not in annotation
            if float(templinelist[headlist.index("hypertrophy_heart_Kin")]) >= rank_thresh:
                return 1
            else:
                return 0
        else:
            return 0
    if phenotype == "failure":
        if templinelist[headlist.index("failure_heart_Kin")] != "": # warning not in annotation
            if float(templinelist[headlist.index("failure_heart_Kin")]) >= rank_thresh:
                return 1
            else:
                return 0
        else:
            return 0

# CURRENTLY UNUSED
#finds variants in topologically central location in network according to user defined criteria
def is_central_all(templine, headlist, rank_thresh, phenotype):
    templinelist = templine.split("\t")
    if phenotype == "normal":
        if templinelist[headlist.index("normal_heart_globalK")] != "":
            if float(templinelist[headlist.index("normal_heart_globalK")]) >= rank_thresh:
                return 1
            else:
                return 0
        else:
            return 0
    if phenotype == "hypertrophy":
        if templinelist[headlist.index("hypertrophy_heart_globalK")] != "":
            if float(templinelist[headlist.index("hypertrophy_heart_globalK")]) >= rank_thresh:
                return 1
            else:
                return 0
        else:
            return 0
    if phenotype == "failure":
        if templinelist[headlist.index("failure_heart_globalK")] != "":
            if float(templinelist[headlist.index("failure_heart_globalK")]) >= rank_thresh:
                return 1
            else:
                return 0
        else:
            return 0

# CURRENTLY UNUSED
#finds variants in expressed regions according to user-defined criteria
def is_expressed(templine, headlist, rank_thresh, phenotype):
    templinelist = templine.split("\t")
    if phenotype == "normal":
        if templinelist[headlist.index("normal_heart_exprs_rank")] != "":
            if float(templinelist[headlist.index("normal_heart_exprs_rank")]) >= rank_thresh:
                return 1
            else:
                return 0
        else:
            return 0
    if phenotype == "hypertrophy":
        if templinelist[headlist.index("hypertrophy_heart_exprs_rank")] != "":
            if float(templinelist[headlist.index("hypertrophy_heart_exprs_rank")]) >= rank_thresh:
                return 1
            else:
                return 0
        else:
            return 0
    if phenotype == "failure":
        if templinelist[headlist.index("failure_heart_exprs_rank")] != "":
            if float(templinelist[headlist.index("failure_heart_exprs_rank")]) >= rank_thresh:
                return 1
            else:
                return 0
        else:
            return 0

# CURRENTLY UNUSED
#finds variants in differentially expressed regions according to user-defined criteria
def is_diffexpr(templine, headlist, phenotype, q_thresh):
    templinelist = templine.split("\t")
    if phenotype == "hypertrophy":
        if templinelist[headlist.index("sam_q_normal_hypertrophy")] != "":
            if float(templinelist[headlist.index("sam_q_normal_hypertrophy")]) <= q_thresh:
                return 1
            else:
                return 0
        else:
            return 0
    if phenotype == "failure":
        if templinelist[headlist.index("sam_q_normal_failure")] != "":
            if float(templinelist[headlist.index("sam_q_normal_failure")]) <= q_thresh:
                return 1
            else:
                return 0
        else:
            return 0

#gives max allele frequency in a list of allele frequency annotations
def max_af(templine, headlist, bp_indexlist):
    templinelist = templine.split("\t")
    af = 0.00
    for i in bp_indexlist:
        if templinelist[i] != "":
            if float(templinelist[i]) > af:
                af = float(templinelist[i])
    return af
                 
#parses info field of vcf file and returns tuple float for mq and mq0
def parse_info(infofield):
    infolist = infofield.split(";")
    for element in infolist:
        if "MQ=" in element:
            mapq=float(element.replace("MQ=", ""))
        elif "MQ0=" in element:
            mapq0=float(element.replace("MQ0=", ""))
    return mapq, mapq0

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
