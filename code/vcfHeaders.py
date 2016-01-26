# list of certain VCF headers from the annotated output file that are used in vcfUtils - NOTE: many of these are now unused. See modules.yml for the headers that are actually being used for tiering, rare variant identification, etc.
# Author: Prag Batra (prag@stanford.edu)

# consts
kCLINVAR_CLNREVSTAT_NO_ASSERTION = 'no_assertion_criteria_provided'
kCLINVAR_CLNREVSTAT_YES_ASSERTION = 'criteria_provided'
kCLINVAR_CLNREVSTAT_ONE_SUBMITTER = 'single_submitter'
kCLINVAR_CLNREVSTAT_CONFLICTING_INTERPS = 'conflicting_interpretations'
# TODO maybe pull the below codes from their definition in the clinvar VCF on import
# new codes for clinvar xml
kCLINVAR_CLINSIG_NOT_PROVIDED_CODE = 'Uncertain significance'
kCLINVAR_CLINSIG_BENIGN_CODE = 'Benign'
kCLINVAR_CLINSIG_LIKELY_BENIGN_CODE = 'Likely benign'
kCLINVAR_CLINSIG_LIKELY_PATHOGENIC_CODE = 'Likely pathogenic'
kCLINVAR_CLINSIG_PATHOGENIC_CODE = 'Pathogenic'
kCLINVAR_CLINSIG_DRUG_RESPONSE_CODE = 'drug response'
kCLINVAR_CLINSIG_OTHER_CODE = 'other'
# old codes for clinvar vcf
# kCLINVAR_CLINSIG_HISTOCOMPATIBILITY_CODE = 7
# kCLINVAR_CLINSIG_NOT_PROVIDED_CODE = 1
# kCLINVAR_CLINSIG_BENIGN_CODE = 2
# kCLINVAR_CLINSIG_LIKELY_BENIGN_CODE = 3
# kCLINVAR_CLINSIG_LIKELY_PATHOGENIC_CODE = 4
# kCLINVAR_CLINSIG_PATHOGENIC_CODE = 5
# kCLINVAR_CLINSIG_DRUG_RESPONSE_CODE = 6
# kCLINVAR_CLINSIG_HISTOCOMPATIBILITY_CODE = 7
# kCLINVAR_CLINSIG_OTHER_CODE = 255

# clinvar headers
kClinvarClinSigHeader = 'clinical_significance'
kClinvarClinRevStatusHeader = 'review_status'
kClinvarStarHeader = 'CLNSTARS' # added by us

# SQL data types for standard VCF columns
# KEYS MUST BE ALL LOWERCASE since we're doing a case insensitive lookup by converting the search term to lowercase
kVCFColTypes = {
                'chrom': 'varchar(10)',
                'start': 'int',
                'pos': 'int', # same as start
                'id': 'varchar(127)',
                'ref': 'varchar(512)',
                'alt': 'varchar(512)',
                'qual': 'varchar(127)',
                'filter': 'varchar(127)'
                # there shouldn't be a definition for the INFO col bc we're only importing specific tags
                }

# for rare variant identification - DEPRECATED and probably UNUSED (now specified in modules YAML instead)
kHapMap2And3_CEU = 'hg19_hapmap2and3_CEU_info'
k1000g_all = 'hg19_popfreq_all_20150413_1000g_all'
k1000g_eur = 'hg19_popfreq_all_20150413_1000g_eur'
kCg69 = 'hg19_cg69_info'
kEsp6500si_ALL = 'hg19_popfreq_all_20150413_esp6500siv2_all'
kEsp6500si_EA = 'hg19_popfreq_all_20150413_esp6500siv2_ea'

k_hapmap2and3_CHB = 'hg19_hapmap2and3_CHB_info'
k_1000g_ASN = '1000g2012apr_ASN' # warning missing
k_hapmap2and3_YRI = 'hg19_hapmap2and3_YRI_info'
k_esp6500si_AA = 'hg19_popfreq_all_20150413_esp6500siv2_aa'
k1000g_afr = 'hg19_popfreq_all_20150413_1000g_afr'

# for tiering - DEPRECATED and probably UNUSED (now specified in modules YAML instead)
kClinvar = 'clinvar_info'
kEncode_dnaseCluster = 'hg19_wgEncodeRegDnaseClustered_r'
kPhyloP_pred = "dbnsfp_PhyloP_Pred" # warning currently not in annotated output file
# kGerp = 'dbnsfp_GERP++' # warning not in annotation
# kGerp = 'GERP++_RS' # check
kPhastConsElts46Way = 'hg19_phastConsElements46way_r_MSA_MCE_score'
# kWg_gerp = 'wg_GERP++_gt2' # warning not in annotation
kWg_gerp = 'hg19_ljb26_all_GERP++_RS'
# kGerp = '' # currently not mapped -- use kWg_gerp instead
kGerp = kWg_gerp # for now
# kSiftPred = 'dbnsfp_SIFT_Pred'
kSiftPred = 'hg19_ljb26_all_SIFT_pred'
# kPolyphen2Pred = 'dbnsfp_PolyPhen2_Pred'
kPolyphen2Pred = 'hg19_ljb26_all_Polyphen2_HVAR_pred' # or polyphen2_HDIV_pred
kPolyphen2Pred_2 = 'hg19_ljb26_all_Polyphen2_HDIV_pred'
# kLrtPred = 'dbnsfp_LRT_Pred'
kLrtPred = 'hg19_ljb26_all_LRT_pred'
# kMutationTasterPred = 'dbnsfp_MutationTaster_Pred'
kMutationTasterPred = 'hg19_ljb26_all_MutationTaster_pred'
