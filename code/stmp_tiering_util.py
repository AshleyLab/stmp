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
import os
import yaml_utils
import yaml_keys
import general_utils
import xlwt # for exporting tiered results as Excel workbook (XLS file)
import stmp_consts
import subprocess

#function for tier 1,2,3,4 targeted calls in specified gene regions, infile is from stanovar vcf
# pop = population our sample comes from
# freq = frequency cutoff
def tiers_allvars(in_file, out_stem, gene_file, pop, yaml_cmds):
    # populate parameters from YAML module specifications
    freq = yaml_cmds[yaml_keys.kModules][yaml_keys.kTiering][yaml_keys.kTRareAlleleFreqCutoff]
    gene_name_col_header = yaml_cmds[yaml_keys.kModules][yaml_keys.kTiering][yaml_keys.kTGeneNameCol]
    functional_column_headers = yaml_utils.convertColumns(yaml_cmds[yaml_keys.kModules][yaml_keys.kTiering][yaml_keys.kTFunctionalCols], yaml_cmds)
    skip_filter_pass_check = yaml_cmds[yaml_keys.kModules][yaml_keys.kTiering][yaml_keys.kTSkipFilterPassCheck]
    
    # check for whether any variant contains "PASS"
    cmd = 'grep "PASS" {infile}|grep -v "#"|wc -l'.format(infile=in_file)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout,stderr = process.communicate()
    returncode = process.returncode
    if(returncode != 0):
        raise ValueError('Failed to search annotated VCF for "PASS" prior to tiering.\ncmd: ' + str(cmd) + '\nErr: ' + str(stderr) + '\nReturncode: ' + returncode + '\nOutput: ' + str(stdout))
    #else
    numHits = int(stdout)
    # TODO move this to pre-checks (prior to even annotation)
    if(numHits == 0 and not skip_filter_pass_check):
        raise ValueError('No variants detected that passed filtering. Re-run STMP with "Skip_Filter_Pass_Check: True" in modules.yml to prioritize all variants anyway.')
    elif(numHits == 0):
        print 'NOTICE: no variants detected that passed filtering. Skipping filter PASS check and prioritizing all variants anyway.'
    else:
        print 'Found ' + str(numHits) + ' variants that passed filtering.' + ' Tiering just these variants.'
        skip_filter_pass_check = False
    
    #open input and output files and initialize counters and lists for background populations
    filein = open(in_file, "r")
    output_log = open(out_stem+".metrics", "w")
    output_log.write("Metrics for stmp filtering, all variants from reference\n")
    header = filein.readline().rstrip("\n")
    headlist = header.split("\t")
    if(gene_file != None):
        g_file = open(gene_file, "r")
    fileoutrare = open(out_stem+'.rare.txt', 'w')
    fileout0 = open(out_stem+".tier0.txt", 'w')
    fileout1 = open(out_stem+".tier1.txt", "w")
    fileout2 = open(out_stem+".tier2.txt", "w")
    fileout3 = open(out_stem+".tier3.txt", "w")
    fileout4 = open(out_stem+".tier4.txt", "w")
    fileoutrare.write(header+"\n")
    fileout0.write("tier\t"+header + "\n")
    fileout1.write("tier\t"+header + "\n")
    fileout2.write("tier\t"+header + "\n")
    fileout3.write("tier\t"+header + "\n")
    fileout4.write("tier\t"+header + "\n")
    total = 0
    damaging0 = 0
    damaging1 = 0
    damaging2 = 0
    damaging3 = 0
    damaging4 = 0
    target_genes = 0
    rarevars = 0
    
    allele_freq_cols = yaml_utils.convertColumns(yaml_cmds[yaml_keys.kModules][yaml_keys.kTiering][yaml_keys.kTAlleleFreqCols], yaml_cmds) #convertTieringColumns(yaml_cmds)
    #debug
    print 'allele freq cols: ' + str(allele_freq_cols)
    
    backpoplist = vcfUtils.get_listindex(headlist, allele_freq_cols)
    #debug
#     print 'backpoplist: ' + str(backpoplist)

    #initialize gene list for region prioritization
    if(gene_file != None):
        genes = {}
        for line in g_file:
            if line.startswith('#'):
                continue
            linelist = line.rstrip("\n").split("\t")
            gene = linelist[0]
            
            if not genes.has_key(gene):
                genes[gene] = 1
            else:
                # debug: uncomment if not debugging
#                 print 'warning: duplicate gene ' + gene + ' in gene list ' + gene_file
                None
    
    #iterate over input file and parse into tiers
    for line in filein:
        total+=1
        if ((skip_filter_pass_check or "PASS" in line) and ("#" in line) == 0 and vcfUtils.is_rare(line, freq, backpoplist, yaml_cmds)
            and not vcfUtils.contains_text('MT', line, [stmp_consts.vcf_col_header_chrom], headlist, yaml_cmds, case_sensitive=False)
            and not vcfUtils.contains_text('ncRNA', line, functional_column_headers, headlist, yaml_cmds, case_sensitive=True)
            ):
            rarevars+=1
            fileoutrare.write(line)
            linelist = line.rstrip("\n").split("\t")
            # for now
            tmp = linelist[headlist.index(gene_name_col_header)].split(',')
            gene = tmp[0]

            if gene_file == None or genes.has_key(gene):
                target_genes+=1
                # tier 0: clinvar
                if(vcfUtils.isClinvarPathogenicOrLikelyPathogenic(line, headlist, yaml_cmds) and not vcfUtils.contains_text('0', line, [yaml_utils.get_datasets(yaml_cmds)['clinvar'][yaml_keys.kDAnnotation]+'_'+vcfHeaders.kClinvarStarHeader], headlist, yaml_cmds, case_sensitive=False)):
                    fileout0.write("0\t"+line)
                    damaging0+=1
                # tier 1
                elif vcfUtils.is_functional(line, "stoploss stopgain splicing frameshift", functional_column_headers, headlist):
                    fileout1.write("1\t"+line)
                    damaging1+=1
                # tier 2
                elif ((vcfUtils.is_functional(line, "nonsynonymous", functional_column_headers, headlist) and vcfUtils.is_conserved(line, headlist, yaml_cmds)) or vcfUtils.is_functional(line, "nonframeshift", functional_column_headers, headlist)):
                    fileout2.write("2\t"+line)
                    damaging2+=1
                # tier 3
                elif vcfUtils.is_functional(line, "nonsynonymous", functional_column_headers, headlist) and vcfUtils.is_pathogenic(line, headlist, yaml_cmds):
                    fileout3.write("3\t"+line)
                    damaging3+=1
                # tier 4
                elif vcfUtils.tolerance_pass(line, headlist, yaml_cmds) and vcfUtils.is_functional(line, "exonic splicing", functional_column_headers, headlist):
                    fileout4.write("4\t"+line)
                    damaging4+=1
                # else ignore variant

    output_log.write("Total variants queried: "+str(total)+"\n")
    output_log.write("Rare variants (allele freq < {freq}) queried: ".format(freq=str(freq))+str(rarevars)+"\n")
    output_log.write("Rare variants in {num} target genes: ".format(num=str(len(genes)) if gene_file != None else '')+str(target_genes)+"\n")
    output_log.write("Candidate variants, tier 0 (rare clinvar pathogenic or likely pathogenic variants with rating > 0 stars): "+str(damaging0)+"\n")
    output_log.write("Candidate variants, tier 1 (rare LOF variants -- stoploss, stopgain, splicing, and frameshift): "+str(damaging1)+"\n")
    output_log.write("Candidate variants, tier 2 (rare nonframeshift or (nonsynonymous and conserved) variants): "+str(damaging2)+"\n")
    output_log.write("Candidate variants, tier 3 (rare nonsynonymous pathogenic variants): "+str(damaging3)+"\n")
    output_log.write("Candidate variants, tier 4 (all other rare exonic/splicing variants with ExAC tolerance z-score (syn_z or mis_z or lof_z) > 2): "+str(damaging4)+"\n")

    filein.close()
    if(gene_file != None):
        g_file.close()
    fileoutrare.close()
    fileout0.close()
    fileout1.close()
    fileout2.close()
    fileout3.close()
    fileout4.close()


# converts tiered output txt files to a single XLS worksheet (with a single sheet per output dir, and variants sorted by tier within each sheet)
def tiers2xls(tier_dirs, output_dir, yaml_commands):
    print 'Saving tiered variants to Excel worksheet'
    # make new excel output worksheet
    wb = xlwt.Workbook(encoding='latin-1') # needed since our SQLite database gives results using latin-1 encoding (so we'd otherwise get UnicodeDecode errors when attempting to decode certain characters if we specified nothing (default is probably ascii) or (probably) UTF-8 as the encoding for the workbook).
    
    for dir in tier_dirs:
        dirname = general_utils.get_file_or_dir_name(dir)
        sheetname = dirname
        sheet = wb.add_sheet(sheetname)
        rowNum = 0 # number of row to write to within excel sheet
        tier0file = os.path.join(dir, 'tiering_allvars.tier0.txt')
        tier1file = os.path.join(dir, 'tiering_allvars.tier1.txt')
        tier2file = os.path.join(dir, 'tiering_allvars.tier2.txt')
        tier3file = os.path.join(dir, 'tiering_allvars.tier3.txt')
        tier4file = os.path.join(dir, 'tiering_allvars.tier4.txt')
        tierFiles = [tier0file, tier1file, tier2file, tier3file, tier4file]
        
        for idx,tierFile in enumerate(tierFiles):
            tierh = open(tierFile, 'r')
            if(idx != 0): # skip header for all files except the first one
                tierh.readline()
                
            warned_of_tier_variants_exceeded = False
            
            for lineNum,line in enumerate(tierh):
                lineContents = line.rstrip("\n").split("\t")
                for col,value in enumerate(lineContents):
                    if(len(value) > 32767):
                        print 'warning: tierFile ' + str(tierFile) + ' line ' + str(lineNum) + ' col ' + str(col) + ' has length ' + str(len(value)) + ' longer than 32767 chars. Truncating.'
                        sheet.write(rowNum, col, value[:32767])
                    elif(lineNum > yaml_commands[yaml_keys.kModules][yaml_keys.kTiering][yaml_keys.kTMaxNumVariantsPerTier]):
                        if(not warned_of_tier_variants_exceeded):
                            print 'warning: number of variants exceeded for ' + str(tierFile) + '. Omitting remaining variants in this tier from XLS file. Current max variants per tier: ' + str(yaml_commands[yaml_keys.kModules][yaml_keys.kTiering][yaml_keys.kTMaxNumVariantsPerTier]) + '. To output more variants per tier, modify this value in modules.yml.'
                            warned_of_tier_variants_exceeded = True
                        continue # don't increment the rowNum since we didn't write anything
                    else:
                        sheet.write(rowNum, col, value)
                
                rowNum += 1
    
    #save the excel file
    out_path = os.path.join(output_dir, 'tiered_output.xls')
    wb.save(out_path)
    return out_path

    
#filter all by local allele frequency
def filter_sfs(infile, sfs_file, outfile, thresh):
    f_sfs = open(sfs_file, "r")
    f_in = open(infile, "r")
    f_out = open(outfile, "w")
    all_dict = {}

    #read in 
    for line in f_sfs:
        linelist = line.replace("\n", "").split("\t")
        chrom = linelist[0].lstrip('chr') # remove chr prefix if present
        #debug
#         print 'sfs chrom: ' + str(chrom)
        pos = linelist[1]
        #debug
#         print 'sfs pos: ' + str(pos)
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
                #debug
#                 print 'tsv chr: ' + str(chrom)
                pos = linelist[posindex]
                #debug
#                 print 'tsv pos: ' + str(pos)
                if all_dict.has_key(chrom+":"+pos):
                    if int(all_dict[chrom+":"+pos]) <= int(thresh):
                        f_out.write(line.replace("\n", "\t")+all_dict[chrom+":"+pos]+"\n")
                else:
                    f_out.write(line.replace("\n", "\t")+"0\n")
    
    f_sfs.close()
    f_in.close()
    f_out.close()
    
