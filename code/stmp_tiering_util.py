#!/usr/bin/env python

'''
Author: Prag Batra prag@stanford.edu

Purpose:
    
    All the heavy lifting for tiering happens here.

Explanation:

    Classifies variant annotations into different tiers of potential pathogenicity according to various criteria.

Example:

    For usage, see stmp.py where these methods are all called.

'''

import vcfUtils
import vcfHeaders
import sys

#function for tier 1,2,3,4 targeted calls in specified gene regions, infile is from stanovar vcf
# pop = population our sample comes from
# freq = frequency cutoff
def tiers_allvars(in_file, out_stem, gene_file, pop, yaml_cmds, freq=0.01, geneNameCol=82):
    #open input and output files and initialize counters and lists for background populations
    filein = open(in_file, "r")
    output_log = open(out_stem+".metrics", "w")
    output_log.write("Metrics for stmp filtering, all variants from reference\n")
    header = filein.readline().rstrip("\n")
    headlist = header.split("\t")
    g_file = open(gene_file, "r")
    fileout0 = open(out_stem+".tier0.txt", 'w')
    fileout1 = open(out_stem+".tier1.txt", "w")
    fileout2 = open(out_stem+".tier2.txt", "w")
    fileout3 = open(out_stem+".tier3.txt", "w")
    fileout4 = open(out_stem+".tier4.txt", "w")
    fileout0.write(header + "\n")
    fileout1.write(header + "\n")
    fileout2.write(header + "\n")
    fileout3.write(header + "\n")
    fileout4.write(header + "\n")
    total = 0
    damaging0 = 0
    damaging1 = 0
    damaging2 = 0
    damaging3 = 0
    damaging4 = 0
    target_genes = 0
    rarevars = 0
    
    if (pop == "CEU") or (pop == "c"):
        backpoplist = vcfUtils.get_listindex(headlist, [vcfHeaders.kHapMap2And3_CEU, vcfHeaders.k1000g_all, vcfHeaders.k1000g_eur, vcfHeaders.kCg69, vcfHeaders.kEsp6500si_ALL, vcfHeaders.kEsp6500si_EA])
    elif (pop == "ASN") or (pop == "a"):
        backpoplist = vcfUtils.get_listindex(headlist, [vcfHeaders.k_hapmap2and3_CHB, vcfHeaders.k1000g_all, vcfHeaders.kCg69, vcfHeaders.kEsp6500si_ALL])
    elif (pop == "AFR") or (pop == "f"):
        backpoplist = vcfUtils.get_listindex(headlist, [vcfHeaders.k_hapmap2and3_YRI, vcfHeaders.k1000g_all, vcfHeaders.k1000g_afr, vcfHeaders.kCg69, vcfHeaders.kEsp6500si_ALL, vcfHeaders.k_esp6500si_AA])
    else:
        print >> sys.stderr, "Error in diseaseUtils.tiers_allvars - Population specified is not supported"
        exit(1)

    #initialize gene list for region prioritization
    genes = {}
    for line in g_file:
        linelist = line.split("\t")
        gene = linelist[1]
        
        if genes.has_key(gene) == 0:
            genes[gene] = linelist[2]+":"+linelist[3]
        else:
            genes[gene] = genes[gene]+";"+linelist[2]+":"+linelist[3]
    
    #iterate over input file and parse into tiers
    for line in filein:
        total+=1
        if ("PASS" in line) and ("#" in line) == 0 and vcfUtils.is_rare(line, freq, backpoplist):
            rarevars+=1
            linelist = line.split("\t")
            # for now
            tmp = linelist[geneNameCol].split(',')
            gene = tmp[0]

            if genes.has_key(gene):
                target_genes+=1
                # tier 0: clinvar
                if(vcfUtils.isClinvarPathogenicOrLikelyPathogenic(line, headlist, yaml_cmds)):
                    fileout0.write(line)
                    damaging0+=1
                elif vcfUtils.is_functional(line, "stoploss stopgain splicing frameshift"):
                    fileout1.write(line)
                    damaging1+=1           
                elif (vcfUtils.is_functional(line, "nonsynonymous") and vcfUtils.is_conserved(line, headlist, 2)) or ("nonframeshift" in line):
                    fileout2.write(line)
                    damaging2+=1
                elif vcfUtils.is_functional(line, "nonsynonymous") and vcfUtils.is_pathogenic(line, headlist, 2):
                    fileout3.write(line)
                    damaging3+=1
                else:
                    fileout4.write(line)
                    damaging4+=1

    output_log.write("Total variants queried: "+str(total)+"\n")
    output_log.write("Rare variants queried: "+str(rarevars)+"\n")
    output_log.write("Rare variants in target genes: "+str(target_genes)+"\n")
    output_log.write("Candidate variants, tier 0: "+str(damaging0)+"\n")
    output_log.write("Candidate variants, tier 1: "+str(damaging1)+"\n")
    output_log.write("Candidate variants, tier 2: "+str(damaging2)+"\n")
    output_log.write("Candidate variants, tier 3: "+str(damaging3)+"\n")
    output_log.write("Candidate variants, tier 4: "+str(damaging4)+"\n")

    filein.close()
    g_file.close()
    fileout0.close()
    fileout1.close()
    fileout2.close()
    fileout3.close()
    fileout4.close()


#function tier 1,2,3,4 variants from reference in specified gene regions, default clinvar - starts from stanovar vcf
def tiers_target(in_file, out_stem, gene_file, pop, yaml_cmds, freq=0.01):
    #open input and output files and initialize counters and lists for background populations
    filein = open(in_file, "r")
    g_file = open(gene_file, "r")
    header = filein.readline().rstrip("\n")
    headlist = header.split("\t")
    output_log = open(out_stem+".metrics", "w")
    output_log.write("Metrics for stmp filtering, clinvar variants\n")
    fileout0 = open(out_stem+".tier0.txt", 'w')
    fileout1 = open(out_stem+".tier1.txt", "w")
    fileout2 = open(out_stem+".tier2.txt", "w")
    fileout3 = open(out_stem+".tier3.txt", "w")
    fileout4 = open(out_stem+".tier4.txt", "w")
    fileout1.write(header + "\n")
    fileout2.write(header + "\n")
    fileout3.write(header + "\n")
    fileout4.write(header + "\n")
    total = 0
    damaging0 = 0
    damaging1 = 0
    damaging2 = 0
    damaging3 = 0
    damaging4 = 0
    target_genes = 0

    # TODO: generate these headers as dynamically as possible using info in the YAML.
    # TODO test to make sure each of these headers actually exists in the annotated file and warn if any are missing
    if (pop == "CEU") or (pop == "c"):
        backpoplist = vcfUtils.get_listindex(headlist, "hg19_hapmap2and3_CEU_info hg19_popfreq_all_20150413_1000g_all hg19_popfreq_all_20150413_1000g_eur hg19_cg69_info hg19_esp6500si_all_info hg19_popfreq_all_20150413_esp6500siv2_all hg19_esp6500si_ea_info hg19_popfreq_all_20150413_esp6500siv2_ea")
#         backpoplist = vcfUtils.get_listindex(headlist, "hapmap2and3_CEU 1000g2010nov_ALL 1000g2011may_ALL 1000g2012apr_ALL 1000g2012apr_EUR cg69 esp6500si_ALL esp6500si_EA")
    elif (pop == "ASN") or (pop == "a"):
        backpoplist = vcfUtils.get_listindex(headlist, "hg19_hapmap2and3_CHB_info hg19_popfreq_all_20150413_1000g_all 1000g2012apr_ASN hg19_cg69_info hg19_esp6500si_all_info hg19_popfreq_all_20150413_esp6500siv2_all") # WARNING missing 1000g_ASN in curent datasets
#         backpoplist = vcfUtils.get_listindex(headlist, "hapmap2and3_CHB 1000g2010nov_ALL 1000g2011may_ALL 1000g2012apr_ALL 1000g2012apr_ASN cg69 esp6500si_ALL")
    elif (pop == "AFR") or (pop == "f"):
        backpoplist = vcfUtils.get_listindex(headlist, "hg19_hapmap2and3_YRI_info hg19_popfreq_all_20150413_1000g_all hg19_popfreq_all_20150413_1000g_afr hg19_cg69_info hg19_esp6500si_all_info hg19_popfreq_all_20150413_esp6500siv2_all hg19_esp6500si_aa_info hg19_popfreq_all_20150413_esp6500siv2_aa")
    else:
        print >> sys.stderr, "Error in diseaseUtils.tiers_allvars - Population specified is not supported"
        exit(1)

    #initialize gene list for region prioritization
    genes = {}
    for line in g_file:
        linelist = line.split("\t")
        gene = linelist[1]
        if genes.has_key(gene) == 0:
            genes[gene] = linelist[2]+":"+linelist[3]
        else:
            genes[gene] = genes[gene]+";"+linelist[2]+":"+linelist[3]

    #iterate over input file and parse into tiers
    for line in filein:
        total+=1
        if ("PASS" in line) and (("#" in line) == 0) and (("0/1" in line) or ("1/1" in line)):
            linelist = line.split("\t")
            gene = linelist[1]
            if genes.has_key(gene):
                target_genes+=1
                # tier 0: clinvar
                if(vcfUtils.isClinvarPathogenicOrLikelyPathogenic(line, headlist, yaml_cmds)):
                    fileout0.write(line)
                    damaging0+=1
                elif vcfUtils.is_functional(line, "stoploss stopgain splicing frameshift"):
                    fileout1.write(line)
                    damaging1+=1
                elif vcfUtils.is_rare(line, freq, backpoplist):
                    fileout2.write(line)
                    damaging2+=1
                elif vcfUtils.is_functional(line, "nonframeshift nonsynonymous"):
                    fileout3.write(line)
                    damaging3+=1
                else:
                    fileout4.write(line)
                    damaging4+=1

    output_log.write("Total variants queried: "+str(total)+"\n")
    output_log.write("Total variants in gene list: "+str(target_genes)+"\n")
    output_log.write("Candidate variants, tier 0: "+str(damaging0)+"\n")
    output_log.write("Candidate variants, tier 1: "+str(damaging1)+"\n")
    output_log.write("Candidate variants, tier 2: "+str(damaging2)+"\n")
    output_log.write("Candidate variants, tier 3: "+str(damaging3)+"\n")
    output_log.write("Candidate variants, tier 4: "+str(damaging4)+"\n")

    filein.close()
    g_file.close()
    fileout0.close()
    fileout1.close()
    fileout2.close()
    fileout3.close()
    fileout4.close()
    
    
#filter all by local allele frequency
def filter_sfs(infile, sfs_file, outfile, thresh):
    f_sfs = open(sfs_file, "r")
    f_in = open(infile, "r")
    f_out = open(outfile, "w")
    all_dict = {}

    #read in 
    for line in f_sfs:
        linelist = line.replace("\n", "").split("\t")
        chrom = linelist[0]
        pos = linelist[1]
    all_dict[chrom+":"+pos] = linelist[2] 

    #read and filter input file
    header = f_in.readline().replace("\n", "\t")+"Num_subjects\n"
    headlist = header.split("\t")
    chrindex = headlist.index("CHROM")
    posindex = headlist.index("POS")
    f_out.write(header)
    for line in f_in:
        if ("CHROM" in line) == 0:
                linelist = line.split("\t")
                chrom = linelist[chrindex]
                pos = linelist[posindex]
                if all_dict.has_key(chrom+":"+pos):
                    if int(all_dict[chrom+":"+pos]) <= int(thresh):
                        f_out.write(line.replace("\n", "\t")+all_dict[chrom+":"+pos]+"\n")
                else:
                    f_out.write(line.replace("\n", "\t")+"0\n")
    
    f_sfs.close()
    f_in.close()
    f_out.close()
    
