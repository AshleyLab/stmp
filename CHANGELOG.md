# Release History

## Version 1.3 (January 8, 2015)
## Highlights
- STMP can now be run in a "test mode" to ensure it is set up correctly. Simply run `python stmp.py --test` to activate this mode. Output will be stored in `(stmp_dir)/data/test/output_unverified` by default, or can be stored in a different directory with the `--output_dir` parameter. If desired, you can compare it against our output (`output_verified`) to see if it is the same (e.g. `diff -rq (stmp_dir)/data/test/output_unverified/ (stmp_dir)/data/test/output_verified/`).
- STMP now includes SnpEFF annotations as part of the joined annotated output file (previously they were only output as a separate file). Note that the SnpEFF annotated column contains the entire VCF INFO column with added SnpEFF annotation tags.
- Order of point and range annotations in the annotated output now reflects the order of the datasets in datasets.yml. Currently, we output VCF columns, Annovar and SnpEFF functional annotations, summary columns, allele frequencies, and other datasets (in that order).
- Added preflight integrity checks to verify that datasets.yml, modules.yml, and the sqlite database file are consistent.
- Tiering now includes all variants only when no variants have "PASS" under FILTER. If some variants have "PASS" in this column, STMP will by default include only these variants in the tiered output. This behavior can be changed by modifyng the Skip_Filter_Pass_Check_If_Needed parameter under tiering in modules.yml.
- Running stmp.py with the --suppress_incidentals flag now causes STMP to ignore clinical gene lists during tiering. As a result, no tabs containing the text "clinical" are present in the tiered output (both folders with text files and Excel spreadsheet).
- The maximum number of variants per tier can now be specified in modules.yml. If any tier contains more than this number of variants, a warning is issued and remaining variants in that tier are suppressed from the Excel output (but are still present in the text output files). By default this number is set to an arbitrarily large value (10,000) to prevent premature truncation of variants. Making this value too large (such that the total number of variants across all tiers exceeds the maximum number of rows allowed in an Excel worksheet) will cause the Excel output to fail.
- Tiering criteria for tier 4 has been tightened to reduce the number of variants in the tiered output. The new tiering criteria is output as part of ".metrics" output file in each tiering folder (e.g. Global, Clinical_Actionable, etc.) and is summarized below for reference (with the change from version 1.2 highlighted in bold):
  - tier 0: rare clinvar pathogenic or likely pathogenic variants
  - tier 1: rare LOF variants -- stoploss, stopgain, splicing, and frameshift
  - tier 2: rare nonframeshift or (nonsynonymous and conserved) variants
  - tier 3: rare nonsynonymous pathogenic variants
  - tier 4: all other rare **exonic/splicing** variants with ExAC tolerance z-score (syn_z or mis_z or lof_z) > 2
- The Import_If_Missing parameter has been added to each dataset in datasets.yml to clarify which datasets should be automatically loaded into the database. Set this parameter to True for all datasets you wish to automatically import on the next database update (stmp.py --update_db), and False for all datasets you do not wish to import at this time. If the False datasets are currently present in the database, their specifications in the database will be checked against datasets.yml to ensure they are consistent.
- Gene description has been added to the default list of consensus (summary) columns in the annotated output. Modifications to the list of summary columns can be made by modifying the Consensus_Columns section of modules.yml.
- Yaml configuration is now exported to the output directory, making it easy to rerun STMP with the same parameters. Note that though the datasets.yml and modules.yml in the output directory are functionally identical to the original configuration files, they do not contain any user-added comments (lines prefixed with "#").
- STMP internal directory structure has been reorganized. This repository now includes stub folders for each subdirectory within STMP, not just the code folder. This should allow for a more streamlined install process (please see the README for updated install instructions).
- stmp.py has been refactored so it can be imported into other Python files without automatically executing the main code.

## Bug fixes and other improvements
- Fixed UnicodeDecode issue which could occur when writing tiered output to Excel file.
- Updated specification for refseq (exonStarts and exonEnds columns are now removed from the database)
- Performance/Debugging: Improved core logic to skip multi-allelic splitting if the file already exists and we are asking for either just the point annotation SQL query or just the range annotation commands.
- Fixed a bug involving datasets ending with "r" (e.g. clinvar).
- Updated SQL query used for point annotations to remove extraneous constraints.
- Improved line count checking for range annotation (so it throws an error if any file has 0 lines).


## Version 1.2 (December 11, 2015)
## Highlights
- Tiered output is now given as an Excel worksheet as well as the previous text files.
- Improved algorithm for functional (Function and ExonicFunction) summary columns. Now the most severe functional annotation is given in each summary column (instead of going down the list of functional datasets and choosing the first dataset that has a valid functional annotation).
- (_experimental_) Added site frequency spectrum (SFS) support.
- Added specifications for gene descriptions and disorder descriptions (for OMIM disease genes) to datasets.yml. (Note: preprocessing required)
- Improved performance of range annotation by separating intersected bed (.tab) file expansion from merging and parallelizing the expansion step.
- Improved tiering algorithm to handle VCFs with all variants missing "PASS" in the FILTER column (by default, we tier all variants instead of just those with "PASS" in the FILTER column. This behavior can be changed in modules.yml under tiering.).
- Updated Bedtools dependency: we are now using version 2.25.0.

## Bug fixes and other improvements
- Fixed several bugs affecting range annotation, particularly for multisample VCFs and VCFs with different numbers of sample columns.
- Removed id, qual, filter columns from the annotated output for VCF datasets.
- Added more information to help debug potential unicode/ascii issues on certain VCFs.
- Improved automatic sanity checking.
- Improved logging of specific steps being performed (for debugging and informational purposes).
- Range annotation: Added check for number of lines in each intersected bed file prior to merging them.
- Annotation: Improved checking of number of lines when merging annotated output files.
- (_experimental_) Wrote sort and merge workaround for range annotation on Mac OS (OSX). This is disabled by default, but can be enabled by uncommenting and commenting the appropriate lines in stmp_annotation_util.py under generate_bed_formatted_dbs.
- Added ability to skip annovar functional annotations (for testing/debugging purposes).


## Version 1.1 (November 17, 2015)
## Highlights
- Annotation and tiering parameters are now stored in modules.yml (located within code/config) for ease of reference and modification. Dataset columns can be specified via either dot notation (yaml dataset name).(column name) or as absolute annotation headers (dataset annotated name)_r(if range annotation)_(column name).
- All tiered output is now stored in a single Excel workbook (tiered_output.xls) in addition to the separate text files described in the README.
- Added support for multimode tiering based on a desired target gene list.By default, runs tiering in 4 different modes: Global (all genes), Clinical_Clinvar (clinvar genes), Clinical_Actionable (acmg+nichorson_114), Clinical_Panel (genedx genes). Output from all of these modes is stored in a single Excel workbook (tiered_output.xls), one sheet per tiering category (e.g. Global, Clinical_Clinvar, etc.), and ordered by tier (tier 0 first). Each tiered output is also stored in a separate subfolder within the output directory. These modes can be modified/removed and additional gene lists can be added in modules.yml (under tiering -> Target_Gene_Lists).
- Added Excel output for tiering (see above). This requires xlwt as a new Python dependency.
- Each row in the tiered output contains the tier number in the first column (e.g. 0, 1, 2, 3, 4).
- Tiering metrics output file now describes the criteria used for each tier as opposed to just giving the tier number.
- *Tiering algorithm changes:* Tier 4 now filters by ExAC tolerance Z-scores to substantially reduce the number of reported variants. All tiers now exclude ncRNA and mitochondrial variants. Tier 0 now only includes clinvar pathogenic/likely pathogenic variants with rating > 0 stars.
- Running stmp.py --update_db now skips previously imported datasets by default (instead of overwriting them). Use --force_overwrite if an existing dataset has changed and all datasets in the YAML will be reimported (NOTE: to save time, it is recommended that you do this with a custom datasets.yml that contains only the datasets you want to overwrite instead of all of them).
- All rare variants (minus ncRNA and mitochondrial variants) are now output in a separate TSV file, in addition to the specific rare variants output as part of tiering.
- Annotated output now contains summary (consensus) columns for gene name, function, and max allele frequency, which show up at the beginning of the annotated output (right after the standard VCF columns). Additional consensus columns can be specified in modules.yml (under annotation -> Consensus_Columns and annotation -> Consensus_Criteria).
- When specifying local dataset files in datasets.yml, files stored in the Dataset_Path (specified in datasets.yml under default) can be referenced by name in addition to absolute path.
- (_experimental_) snpEFF annotations are output as a separate TSV file (scratch.genes.txt).
- Swapped original Clinvar VCF for MacArthur VCF generated from parsing Clinvar XML as the original VCF had multiple inaccuracies.

## Bug fixes and other improvements
- By default, tiering no longer filters by gene list (previously we were filtering by an outdated ClinVar gene list).
- YAML parsing code and constants, as well as general have been moved to separate python files for greater modularity.
- Additional code refactoring for greater modularity.
- Fixed issue involving case-sensitivity in imported dataset/column names and another issue involving modified names for region datasets.
- Improved tiering checks to reduce potential false positives.
- Logging: Added basic checks on the annotated output file to warn about any columns that are empty for all variants.
- Path for BED files (used in range annotation) is now specified in datasets.yml (under default -> Bed_Path).
- BED files for range annotation are now generated on db update instead of the first annotation run to avoid potential issues with multiple concurrent annotation runs immediately after db update.
- Dataset downloads are now performed using wget --no-check-certificate, which forces the file download even if the website certificates don't match (e.g. for github file downloads).
- Improved clinvar star annotation (and made it work with the new MacArthur Clinvar VCF above).
- Fixed bug involving IndexError on range annotation with larger input files.
- Improved logging and error reporting if a numeric column used in tiering (e.g. allele frequency) contains non-numeric data.
- Improved logging for third party tools (added print statements for range annotation and other third party tool invocations).
