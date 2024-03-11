# code_release_001

In here you will find the code that was used to analyse the species information coming from the analysis of the ARMS-MBON sequences

* [WormsAttributes4ARMSdata.py](https://github.com/arms-mbon/code_release_001/blob/main/WormsAttributes4ARMSdata.py) takes as input a CSV file with at least one colunm of [AphiaIDs](https://www.marinespecies.org/about.php#what_is_aphia) for the species' of interest, and it returns the information about a set of attributes ("Species importance to society", "IUCN RedList Category","IUCN Criteria","IUCN Year Accessed","HELCOM RedList Category","AMBI ecological group","Environmental position") for those species as found in WoRMS, using its [REST APIs](https://www.marinespecies.org/rest/). Useage of the code is documented within the code. You run this on the command line, with the input file name and column number with the AphiaIDs written into the top of the code.   
