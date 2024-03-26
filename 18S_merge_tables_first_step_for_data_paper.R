library(tidyr)
library(tibble)
library(Biostrings)
library(dplyr)                                                 
library(plyr)                                                  
library(readxl)
library(stringr)
library(openxlsx)    

setwd("~/18S")

# read Extended_final_tables of each sequencing run
# One Gdynia sample was processed separately and therefore has its own table

OTUcounts_April2021 <- read_xlsx("Extended_final_table_April2021_18S_noBlank.xlsx")

OTUcounts_Jan2020 <- read_xlsx("Extended_final_table_January2020_18S_noBlank.xlsx")

OTUcounts_Jan2022 <- read_xlsx("Extended_final_table_January2022_18S_noBlank.xlsx")

OTUcounts_July2019 <- read_xlsx("Extended_final_table_July2019_18S_noBlank.xlsx")

OTUcounts_May2021 <- read_xlsx("Extended_final_table_May2021_18S_noBlank.xlsx")

OTUcounts_Sep2020 <- read_xlsx("Extended_final_table_September2020_18S_noBlank.xlsx")

OTUcounts_Aug2023 <- read_xlsx("Extended_final_table_August2023_18S_noBlank.xlsx")

OTUcounts_GDY <- read_xlsx("Extended_final_table_ARMS_Gdynia_GDY1_20180813_20191029_SF38_DMSO_18S_noBlank.xlsx")

# read fasta files to get sequences 
# all_sequences_grouped_XX_18S.fa files contain all clustered OTU sequences, including OTUs present in blank samples and OTUs not classified taxonomically
# the following fasta files therefore contain more sequences than OTUs present in the Extended_final_table_XX_18S_noBlank.xlsx files

fasta_April2021 <- readDNAStringSet("all_sequences_grouped_April2021_18S.fa")

fasta_Jan2020 <- readDNAStringSet("all_sequences_grouped_January2020_18S.fa")

fasta_Jan2022 <- readDNAStringSet("all_sequences_grouped_January2022_18S.fa")

fasta_July2019 <- readDNAStringSet("all_sequences_grouped_July2019_18S.fa")

fasta_May2021 <- readDNAStringSet("all_sequences_grouped_May2021_18S.fa")

fasta_Sep2020 <- readDNAStringSet("all_sequences_grouped_September2020_18S.fa")

fasta_Aug2023 <- readDNAStringSet("all_sequences_grouped_August2023_18S.fa")

fasta_GDY <- readDNAStringSet("all_sequences_grouped_ARMS_Gdynia_GDY1_20180813_20191029_SF38_DMSO_18S.fa")

# Make lists of OTUcounts and fastas

all_OTUcounts<-tibble::lst(OTUcounts_April2021,OTUcounts_Jan2020,OTUcounts_Jan2022,OTUcounts_July2019,OTUcounts_May2021,OTUcounts_Sep2020,OTUcounts_Aug2023,OTUcounts_GDY) # lst from tibble will preserve object names
all_fasta<-tibble::lst(fasta_April2021,fasta_Jan2020,fasta_Jan2022,fasta_July2019,fasta_May2021,fasta_Sep2020,fasta_Aug2023,fasta_GDY) 

# Make fasta files into dataframes 

IDs<-list()
sequences<-list()
all_seqs<-list()
for(i in 1:length(all_fasta)) {
  IDs[[i]]<-names(all_fasta[[i]])
  sequences[[i]]<-paste(all_fasta[[i]])
  all_seqs[[i]]<-data.frame(IDs[[i]],sequences[[i]])
  colnames(all_seqs[[i]])<-c("ID","seqs")
}  

# Get number of unique OTUs prior to any curation

length(unique(rbind(all_seqs[[1]],all_seqs[[2]],all_seqs[[3]],all_seqs[[4]],all_seqs[[5]],all_seqs[[6]],all_seqs[[7]],all_seqs[[8]])$seqs))

# Subset to OTUs found in count tables (these are the OTUs with a taxonomic classification remaining after blank-curation)

for(i in 1:length(all_seqs)) {
  all_seqs[[i]] <- all_seqs[[i]] %>% filter(ID %in% all_OTUcounts[[i]]$OTU)
}  

# Add sequences to OTUcounts

for(i in 1:length(all_OTUcounts)) {
  for(j in 1:nrow(all_OTUcounts[[i]])) {
    all_OTUcounts[[i]]<-data.frame(all_OTUcounts[[i]])
    all_OTUcounts[[i]]$seqs[all_OTUcounts[[i]]$OTU == all_seqs[[i]]$ID[j]] <- all_seqs[[i]][,2][j]
  }
}

# Make taxonomy tables and subset to OTUs present in each table of all_OTUcounts 

taxa<-lapply(all_OTUcounts,function(x) x %>% select(c(seqs,Classification)))

# Bring sequence columns to front, remove OTU and taxonomy columns

all_OTUcounts<-lapply(all_OTUcounts,function(x) x %>% relocate(seqs))
all_OTUcounts<-lapply(all_OTUcounts,function(x) x %>% select(-c(OTU,Classification,TAXON.NCBI_TAX_ID)))

# Combine OTUcounts, set NAs to zero

OTUcounts<-rbind.fill(all_OTUcounts[[1]],all_OTUcounts[[2]],all_OTUcounts[[3]],all_OTUcounts[[4]],all_OTUcounts[[5]],all_OTUcounts[[6]],all_OTUcounts[[7]],all_OTUcounts[[8]])
OTUcounts[is.na(OTUcounts)]<-0

# Aggregate OTU counts based on unique sequences

OTUcounts<-aggregate(.~ seqs,data = OTUcounts,FUN=sum)

# Combine taxonomy tables

Taxfile<-rbind(taxa[[1]],taxa[[2]],taxa[[3]],taxa[[4]],taxa[[5]],taxa[[6]],taxa[[7]],taxa[[8]])

# Bring taxonomy strings into separate columns and set certain entries as NA

Taxfile<-separate_wider_delim(Taxfile,Classification, delim = ";", names_sep="",too_few = "align_start")

Taxfile[] <- lapply(Taxfile, function(x) replace(x, grepl("var\\.", x), "var.")) # keep strings containing "var." just as "var."

Taxfile[] <- lapply(Taxfile, function(x) replace(x, grepl(" |X |XX|sp\\.", x), NA)) # sets all levels with a space, certain capital X strings or sp. in the string as NA

Taxfile<-Taxfile[colSums(!is.na(Taxfile)) > 0] # Remove columns which are now left with NAs only

# Keep only distinct sequence-taxonomy combinations

Taxfile<-Taxfile %>% distinct()

# Paste taxonomy levels together ignoring NA entries

Taxfile<-transform(Taxfile, Taxonomy=sapply(apply(Taxfile[,2:ncol(Taxfile)], 1, \(x) x[!is.na(x)]), paste, collapse='_'))
Taxfile<-Taxfile %>% select(c(seqs,Taxonomy))

# Separate taxonomy strings again

Taxfile<-separate_wider_delim(Taxfile,Taxonomy, delim = "_", names_sep="",too_few = "align_start")

# Create a Species column combining species and genus level in one string
# This is done by checking for strings which do not have a capital letter (i..e, species level names)
# Done up until 3rd column just to be sure

Taxfile<-as.data.frame(Taxfile)
Taxfile[Taxfile=="clade"]<-"Clade" # change lowercase clade to uppercase Clade
ncol(Taxfile) # check number of columns
Taxfile$Species<-NA # Make Species column and fill with NA as dummy
Taxfile$Species<-ifelse(str_detect(Taxfile[,10],"[[:upper:]]",negate=T) & !is.na(Taxfile[,10]),paste(Taxfile[,9],Taxfile[,10],sep="_"),Taxfile$Species)
Taxfile$Species<-ifelse(str_detect(Taxfile[,9],"[[:upper:]]",negate=T) & !is.na(Taxfile[,9]),paste(Taxfile[,8],Taxfile[,9],sep="_"),Taxfile$Species)
Taxfile$Species<-ifelse(str_detect(Taxfile[,8],"[[:upper:]]",negate=T) & !is.na(Taxfile[,8]),paste(Taxfile[,7],Taxfile[,8],sep="_"),Taxfile$Species)
Taxfile$Species<-ifelse(str_detect(Taxfile[,7],"[[:upper:]]",negate=T) & !is.na(Taxfile[,7]),paste(Taxfile[,6],Taxfile[,7],sep="_"),Taxfile$Species)
Taxfile$Species<-ifelse(str_detect(Taxfile[,6],"[[:upper:]]",negate=T) & !is.na(Taxfile[,6]),paste(Taxfile[,5],Taxfile[,6],sep="_"),Taxfile$Species)
Taxfile$Species<-ifelse(str_detect(Taxfile[,5],"[[:upper:]]",negate=T) & !is.na(Taxfile[,5]),paste(Taxfile[,4],Taxfile[,5],sep="_"),Taxfile$Species)
Taxfile$Species<-ifelse(str_detect(Taxfile[,4],"[[:upper:]]",negate=T) & !is.na(Taxfile[,4]),paste(Taxfile[,3],Taxfile[,4],sep="_"),Taxfile$Species)
Taxfile$Species<-ifelse(str_detect(Taxfile[,3],"[[:upper:]]",negate=T) & !is.na(Taxfile[,3]),paste(Taxfile[,2],Taxfile[,3],sep="_"),Taxfile$Species)

# Sort Taxfile based on order in OTUcounts

Taxfile<-Taxfile[order(match(Taxfile[,1],OTUcounts[,1])),]

# Set names of taxonomic levels

colnames(Taxfile)[2:(ncol(Taxfile)-1)]<-c("Domain","Supergroup","Division/Kingdom","Phylum/Class","Level_X","Level_XX","Level_XXX","Level_XXXX","Level_XXXXX")

# Set assignments containing "lineage" string in Species column as NA

Taxfile$Species <- replace(Taxfile$Species, grepl("lineage", Taxfile$Species),NA)

# Some cells in Taxfile may have a "var." string remaining. Check if this is the case and manually adjust taxonomy for those cases.

Taxfile_var <- Taxfile[apply(Taxfile, 1, function(x) any(grepl("var.", x))), ]
Taxfile_var # One row with "var." string
Taxfile[Taxfile$seqs==Taxfile_var$seqs,"Level_XXXX"]<-"var. lyococcos"
Taxfile[Taxfile$seqs==Taxfile_var$seqs,"Level_XXXXX"]<-NA
Taxfile[Taxfile$seqs==Taxfile_var$seqs,"Species"]<-"Rhizopus_stolonifer_var._lyococcos"

# remove numeric characters and "-" remaining in species names

Taxfile$Species <- gsub('[0-9.]', '', Taxfile$Species)
Taxfile$Species <- gsub('-', '', Taxfile$Species)

# Make fasta file of unique sequences with short OTU names and save to file

fasta<-data.frame(OTU=paste0(">OTU",1:nrow(OTUcounts)),seqs=OTUcounts$seqs) |>
  pivot_longer(everything()) |> 
  subset(select=-name)

write.table(fasta,"18S_OTUs.fa",row.names = F,col.names = F,quote=F)

# Write OTU table to file

OTUcounts$OTU<-paste0("OTU",1:nrow(OTUcounts))
OTUcounts<-OTUcounts %>% relocate(OTU) %>% select(-seqs)
write.table(OTUcounts,"18S_OTU_table.txt",sep="\t",row.names = F)

# Write taxonomy table to file

Taxfile$OTU<-paste0("OTU",1:nrow(Taxfile))
Taxfile<-Taxfile %>% relocate(OTU) %>% select(-seqs)
write.table(Taxfile,"18S_tax_table.txt",sep="\t",row.names = F)

# get number of unique, classified, blank-corrected OTUs

length(OTUcounts$OTU)

# Get number of reads

sum(colSums(OTUcounts[,-1]))

