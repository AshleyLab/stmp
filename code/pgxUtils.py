#Rick Dewey 2.13.13 rdewey@stanford.edu, Prag Batra 2015 prag@stanford.edu
#PGX annotation
#utils for PGX annotations for STMP

import os
import re
import fileinput
import sys
import getopt
import glob
import gzip


def open_compressed_or_regular(f, options):
    if(f.endswith('.gz')):
        return gzip.open(f, options)
    #else
    return open(f, options)

#reads in vcf file and return variant dictionary
def vcf_reader_vars(infile):
    var_dict = {}
    
    filein = open_compressed_or_regular(infile, "r")
    for line in filein.readlines():
        if ("#" in line) == 0:
            linelist = line.replace("\n", "").split("\t")
            rsid = linelist[2]
            ref = linelist[3]
            alt = linelist[4]
            gt_list = linelist[9].split(":")
            gt = gt_list[0]
            if gt == "0/0":
                var_dict[rsid] = ref+ref
            elif gt == "1/1":
                    var_dict[rsid] = alt+alt
            elif gt == "0/1":
                if ref<alt:
                        var_dict[rsid] = ref+alt
                else:
                    var_dict[rsid] = alt+ref
    return var_dict


#reads in vcf file and return variant dictionary, list of heterozygous variants, list of homozygous variants
def vcf_reader_haps(infile):
    var_dict = {}
    hom_list = []
    het_list = []
    all_list = []
    
    filein = open_compressed_or_regular(infile, "r")
    for line in filein.readlines():
        if ("#" in line) == 0:
            linelist = line.replace("\n", "").split("\t")
            rsid = linelist[2]
            ref = linelist[3]
            alt = linelist[4]
            gt_list = linelist[9].split(":")
            gt = gt_list[0]
            #debug
#             print 'gt: ' + str(gt) + ' rsid: ' + str(rsid) + ' ref: ' + str(ref) + ' alt ' + str(alt) + ' gt_list ' + str(gt_list)
            if gt == "0/0":
                var_dict[rsid] = ref
                hom_list.append(rsid)
                all_list.append(rsid)
            elif gt == "0/1":
                var_dict[rsid] = ref+alt
                het_list.append(rsid)
                all_list.append(rsid)
            elif gt == "1/1":
                var_dict[rsid] = alt
                hom_list.append(rsid)
                all_list.append(rsid)
    return var_dict, hom_list, het_list, all_list


#annotate variants according to pharmgx clinical annotations
def pgx_annotator(n_file, p_file, o_stem):
    pgx_in = open_compressed_or_regular(p_file, "r")
    fout_tox = open_compressed_or_regular(o_stem+".toxicity.txt", "w")
    fout_eff = open_compressed_or_regular(o_stem+".efficacy.txt", "w")
    fout_dos = open_compressed_or_regular(o_stem+".dosage.txt", "w")
    fout_other = open_compressed_or_regular(o_stem+".other.txt", "w")
    pgx_log = open_compressed_or_regular(n_file.replace(".vcf", ".pgx_annotations.log"), "w")
    pgx_log.write("Pharmacogenomics annotator log")
    fout_tox.write("rsid\t"+"gene\t"+"genotype\t"+"annotation\t"+"level of evidence\t"+"annotation type\t"+"PMIDs\t"+"drugs\t"+"disease\t"+"ethnicity\n")
    fout_eff.write("rsid\t"+"gene\t"+"genotype\t"+"annotation\t"+"level of evidence\t"+"annotation type\t"+"PMIDs\t"+"drugs\t"+"disease\t"+"ethnicity\n")
    fout_dos.write("rsid\t"+"gene\t"+"genotype\t"+"annotation\t"+"level of evidence\t"+"annotation type\t"+"PMIDs\t"+"drugs\t"+"disease\t"+"ethnicity\n")
    fout_other.write("rsid\t"+"gene\t"+"genotype\t"+"annotation\t"+"level of evidence\t"+"annotation type\t"+"PMIDs\t"+"drugs\t"+"disease\t"+"ethnicity\n")
    trash = pgx_in.readline()

    efficacy = 0
    dosage = 0
    toxicity = 0
    other = 0

    vcf_vars = vcf_reader_vars(n_file)

    for line in pgx_in.readlines():
        line_list = line.replace('\n', '').replace('"', '').split("\t")
        rsid = line_list[1]
        annot = line_list[6].split(";")
        if vcf_vars.has_key(rsid):
            for item in annot:
                item_list = item.split(":")
                gt = item_list[0].replace('"', '').replace(' ', '')

                #permissive definitions for parsing variant output
                if vcf_vars[rsid] == gt:
                    if "Efficacy" in line_list[4]:
                        fout_eff.write(rsid+"\t"+line_list[2]+"\t"+vcf_vars[rsid]+"\t"+"".join(item_list[1:len(item_list)])+"\t"+line_list[3]+"\t"+line_list[4]+"\t"+line_list[9]+"\t"+line_list[11]+"\t"+line_list[12]+"\t"+line_list[13]+"\n")
                        efficacy+=1
                    if "Dosage" in line_list[4]:
                        fout_dos.write(rsid+"\t"+line_list[2]+"\t"+vcf_vars[rsid]+"\t"+"".join(item_list[1:len(item_list)])+"\t"+line_list[3]+"\t"+line_list[4]+"\t"+line_list[9]+"\t"+line_list[11]+"\t"+line_list[12]+"\t"+line_list[13]+"\n")
                        dosage+=1
                    if "Toxicity" in line_list[4]: 
                        fout_tox.write(rsid+"\t"+line_list[2]+"\t"+vcf_vars[rsid]+"\t"+"".join(item_list[1:len(item_list)])+"\t"+line_list[3]+"\t"+line_list[4]+"\t"+line_list[9]+"\t"+line_list[11]+"\t"+line_list[12]+"\t"+line_list[13]+"\n")
                        toxicity+=1
                    if "Other" in line_list[4]: 
                        fout_other.write(rsid+"\t"+line_list[2]+"\t"+vcf_vars[rsid]+"\t"+"".join(item_list[1:len(item_list)])+"\t"+line_list[3]+"\t"+line_list[4]+"\t"+line_list[9]+"\t"+line_list[11]+"\t"+line_list[12]+"\t"+line_list[13]+"\n")
                        other+=1
        else:
            pgx_log.write("\nrsid not found in vcf file: "+rsid)

    pgx_log.write("\n\nDrug efficacy annotations: "+str(efficacy))
    pgx_log.write("\nDrug dosage annotations: "+str(dosage))
    pgx_log.write("\nDrug toxicity annotations: "+str(toxicity))
    pgx_log.write("\nOther drug response annotations: "+str(other))
    fout_tox.close()
    fout_eff.close()
    fout_dos.close()
    fout_other.close()
    pgx_in.close()
    pgx_log.close()


#for a given set of inputs from vcf_reader and haplotype table, returns a dict with possible star allele combos
def haplotyper(gene, hapfile, vcf_list, vhet_list, vhom_list, vdict, in_file):

    #debug
#     print 'VCF list: ' + str(vcf_list)
#     print 'hapfile: ' + str(hapfile)
#     print 'gene: ' + str(gene)
#     print 'vhet_list: ' + str(vhet_list)
#     print 'vhom_list: ' + str(vhom_list)
#     print 'vdict: ' + str(vdict)
#     print 'infile: ' + str(in_file)

    hap_dict = {}
    alleles = {}
    use_index = []
    rs_list = []
    haplotypes = []
#     logfilename = in_file.replace(".vcf", "_"+gene+".log")
    #debug
#     print 'pgx logfile: ' + logfilename
#     logfile = open_compressed_or_regular(logfilename, "w")
#     logfile = open(sys.stdout, 'w')

    count_het = 0
    
    print "Pharmacogenomics haplotyper log for "+gene
#     logfile.write("Pharmacogenomics haplotyper log for "+gene)

    #read in haplotype file
    hapin = open_compressed_or_regular(hapfile, "r")
    header = hapin.readline()
    head_list = header.replace("\n", "").split("\t")
    var_list = head_list[1:len(head_list)]

    var_set = set(var_list)
    vcf_set = set(vcf_list)
    use_var = var_set.intersection(vcf_set)

    #perform set intersection between confidently called variants in vcf file and variants represented on pharmgkb haplotypes
    for element in list(use_var):
        use_index.append(head_list.index(element))
    not_found = var_set.difference(vcf_set)

    use_index.sort()
    
    if(len(use_index) == 0): # no intersections bw vcf variants and haplotype table
        print 'warning: no intersections found between variants in VCF and hapfile ' + str(hapfile) + ' for gene ' + gene
        alleles = {}
        found_pair = 0
        return alleles, found_pair

    #ordered list of haplotype rsids corresponding to haplotype designation from pharmgkb
    for element in use_index:
        rs_list.append(head_list[element])

    #initialize list for haplotypes, there are 2^n possible haplotypes for n heterozygous variants, and 2^(n-1) possible haplotype combinations
    for rsid in rs_list:
        if (rsid in vhet_list):
            count_het+=1

    print("\n\nNumber of genotypes used for allele approximation: "+str(len(use_index)))
    print("\nNumber of heterozygous variants at pgx loci: "+str(count_het))
    print("\nPgx variants not found in vcf call file: "+",".join(not_found))
    print("\nGenotype calls:")

    if rs_list[0] in vhom_list:
        var_haps = [vdict[rs_list[0]]]
        print("\n"+rs_list[0]+":"+vdict[rs_list[0]]+vdict[rs_list[0]])
    elif rs_list[0] in vhet_list:
        allele1 = vdict[rs_list[0]][0]
        allele2 = vdict[rs_list[0]][1]
        var_haps = [allele1, allele2]
        print("\n"+rs_list[0]+":"+allele1+allele2)

    for element in rs_list[1:len(rs_list)]:
        haps_old = var_haps
        if element in vhom_list:
            var_haps = []
            for i in range(0,len(haps_old)):
                var_haps.append(haps_old[i]+vdict[element])
            print("\n"+element+":"+vdict[element]+vdict[element])
        elif element in vhet_list:
            var_haps = []
            allele1 = vdict[element][0]
            allele2 = vdict[element][1]
            for j in range(0,len(haps_old)):
                var_haps.append(haps_old[j]+allele1)
                var_haps.append(haps_old[j]+allele2)
            print("\n"+element+":"+allele1+allele2)

    hap_pair = {}

    #duplicate haplotype for uniformly homozygous calls
    if len(var_haps) == 1:
        var_haps.append(var_haps[0])

    #create hapolotype pairs
    for j in range(0, len(var_haps)/2):
        hap_pair[var_haps[j]] = var_haps[len(var_haps)-j-1]
    
    #create dictionary of haplotypes from pharmgkb and associated star allele designations
    for line in hapin.readlines():
        haplist = []
        line_list = line.replace("\n", "").split("\t")
        for element in use_index:
            haplist.append(line_list[element])
        haplotype = "".join(haplist)
        if hap_dict.has_key(haplotype):
            hap_dict[haplotype] = hap_dict[haplotype]+","+line_list[0]
        else:
            hap_dict[haplotype] = line_list[0]

    alleles = {}
    found_pair = 0

    #search for haplotype pairs in allele lookup table
    for key, value in hap_pair.iteritems():
        if hap_dict.has_key(key):
            allele1 = key+":"+hap_dict[key]
            if hap_dict.has_key(value):
                allele2 = value+":"+hap_dict[value]
                alleles[allele1] = allele2
                found_pair+=1

    print("\n\nNumber of possible pgx haplotype pairs found: "+str(found_pair))
    if found_pair >0:
        print("\n"+"Possible allele pairs: ")
        for key, value in alleles.iteritems():
            print("\n"+key+","+value)

    return alleles, found_pair


#wrapper for calling star alleles for all haplotype tables in a directory
def star_caller(hap_dir, n_file, o_stem):
    
    #debug
#     print 'hap_Dir: ' + str(hap_dir)
#     print 'n_file: ' + str(n_file)
#     print 'o_stem: ' + str(o_stem)
     
    hap_file_list = glob.glob(hap_dir+"/*-haplotypes-snvs.txt")
    fileout = open_compressed_or_regular(o_stem+".star.txt", "w")

    #read in input vcf file
    vcf_dict, vcf_hom, vcf_het, vcf_all = vcf_reader_haps(n_file)

    #iterate over haplotypes in directory of snv-defined star alleles
    for item in hap_file_list:
        gene_symb = item.replace(hap_dir, "").replace("/", "").split("-")[0]
        all_dict, is_found =  haplotyper(gene_symb, item, vcf_all, vcf_het, vcf_hom, vcf_dict, n_file)
        if is_found > 0:
            for key, value in all_dict.iteritems():
                fileout.write(gene_symb+"\t"+key+"\t"+value+"\n")
        else:
            fileout.write(gene_symb+"\t"+"No haplotype matches\n")
    fileout.close()

