# code_release_001

Here you will find all code that was used for the analysis presented in ARMS-MBON's first data paper. Note that this is not the code used to create files for EurOBIS submissions, but the code used for the exploration of the sequencing data of data_release_001.

## Processing in R

The following R scripts merge the data (read count, taxonomy, and fasta files) from individual PEMA runs we provide in the [analysis_release_001](https://github.com/arms-mbon/analysis_release_001/tree/main) GitHub repository for each marker gene: 

* [COI_merge_tables_first_step_for_data_paper.R](https://github.com/arms-mbon/code_release_001/blob/main/COI_merge_tables_first_step_for_data_paper.R)
* [18S_merge_tables_first_step_for_data_paper.R](https://github.com/arms-mbon/code_release_001/blob/main/18S_merge_tables_first_step_for_data_paper.R)
* [ITS_merge_tables_first_step_for_data_paper.R](https://github.com/arms-mbon/code_release_001/blob/main/ITS_merge_tables_first_step_for_data_paper.R) 

These scripts perform the following tasks:
* As no confidence threshold was applied within PEMA for taxonomic assignments of COI ASVs, all rank assignments with a confidence value of below 0.8 are discarded for this marker gene.
* Given that reference databases use different taxonomic levels and contain assignments that may not represent an actual classification (e.g., "Class.xy_X" etc.), rank assignments were curated (i.e., assignments containing Xs,"_sp", etc. set as NA), actual species assignments generated and a final rank order established.
* Final ASV/OTU count tables, taxonomy tables and fasta files were generated for each gene's data set. These files can be found [here](https://github.com/arms-mbon/code_release_001/tree/main/final_count_taxonomy_fasta_files).

The R script below performs all subsequent exploration of the sequencing data, including 
To analyse the species information coming from the analysis of the ARMS-MBON sequences:
* [WormsAttributes4ARMSdata.py](https://github.com/arms-mbon/code_release_001/blob/main/WormsAttributes4ARMSdata.py) takes as input a CSV file with at least one colunm of [AphiaIDs](https://www.marinespecies.org/about.php#what_is_aphia) for the species' of interest, and it returns the information about a set of attributes ("Species importance to society", "IUCN RedList Category","IUCN Criteria","IUCN Year Accessed","HELCOM RedList Category","AMBI ecological group","Environmental position") for those species as found in WoRMS, using its [REST APIs](https://www.marinespecies.org/rest/). Useage of the code is documented within the code. You run this on the command line, with the input file name and column number with the AphiaIDs written into the top of the code.   
