# constant keys in the yaml file (for config information)
# Author: Prag Batra (prag@stanford.edu)


############## DATASETS ###############

kDatasets = 'Datasets' # key to access unordered dictionary of datasets (with all info) in parsed yaml
kDatasetDefaults = 'Dataset_Defaults' # key to access default parameters in parsed yaml
kDatasetOrder = 'Dataset_Order' # key to access ordered list of dataset NAMES ONLY in parsed yaml

# key for default parameters in datasets.yml (prior to parsing)
kDDefaults = 'default'
kDBuild = 'Build'
kDDatasetsPath = 'Dataset_Path'
# Dataset description keys
kDImportIfMissing = 'Import_If_Missing'
kDBedPath = 'Bed_Path'
kDAnnotation = 'Annotation'
kDCategory = 'Category'
kDCategoryTypeRegion = 'region'
kDMultimatchDelimiter = 'Delimiter_for_multiple_matches'
kDColumnHeaders = 'ColumnHeaders'
kDColumnTypes = 'DataType'
kDSource = 'Source'
kDStartingLine = 'StartingLine'


############## MODULES ###############

kModules = 'Modules'

## Testing keys
kTesting = 'testing'
kTeTestDatasetPath = 'Test_Dataset_Path'
kTeTestDefaultOutPath = 'Default_Test_Output_Dir'

## Annotation keys
kAnnotation = 'annotation' # annotation module spec, NOT dataset parameter
kAAnnovarPath = 'Annovar_Path'
kASnpeffPath = 'Snpeff_Path'
kASnpeffMemory = 'Snpeff_Memory'
kAConsensusColumnSuffix = 'Consensus_Column_Suffix'
kAConsensusColumnOffset = 'Consensus_Column_Insertion_Offset'
kAConsensusColumns = 'Consensus_Columns'
kAConsensusCriteria = 'Consensus_Criteria'
kAConsensusCriteriaTypeMaxFreq = 'Max_Freq'
kAConsensusColumnsOrder = 'Order'
kAConsensusColumnsGene = 'Gene'
kABedMultimatchInternalDelimiter = 'Bed_Multimatch_Internal_Delimiter'
kABedInternalDelimiter = 'Bed_Internal_Delimiter'

## Tiering
kTiering = 'tiering'
kTMaxNumVariantsPerTier = 'Max_Num_Variants_Per_Tier'
kTSkipFilterPassCheck = 'Skip_Filter_Pass_Check_If_Needed'
kTRareAlleleFreqCutoff='Rare_Allele_Frequency_Cutoff'
kTGeneNameCol='Gene_Name_Column'
kTAlleleFreqCols='Allele_Frequency_Columns'
kTFunctionalCols = 'Functional_Columns'
kTColMultipleThresholdSeparator = 'Column_Multiple_Threshold_Separator'
# conservation
kTConservationCols = 'Conservation_Columns'
kTConservationThresholds = 'Conservation_Thresholds'
kTConservationGlobalThreshold = 'Conservation_Global_Threshold'
# pathogenicity
kTPathogenicityCols = 'Pathogenicity_Columns'
kTPathogenicityThresholds = 'Pathogenicity_Thresholds'
kTPathogenicityGlobalThreshold = 'Pathogenicity_Global_Threshold'
# tolerance
kTToleranceZScoreCols = 'Tolerance_Zscore_Columns'
kTToleranceZScoreCutoff = 'Tolerance_Zscore_Cutoff'
# Targeted gene lists
kTTargetGeneLists = 'Target_Gene_Lists'
kTClinicalGeneList = 'Clinical'
kTClinicalClinvarGeneList = 'Clinical_Clinvar'
kTClinicalPanelGeneList = 'Clinical_Panel'
kTActionableGeneList = 'Actionable'


## PGX
# TBA


## Trio
# TBA

