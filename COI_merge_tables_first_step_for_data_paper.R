library(tidyr)
library(tibble)
library(Biostrings)
library(dplyr)                                                 
library(plyr)                                                  
library(readxl)
library(openxlsx)

setwd("C:/Users/Nauras/OneDrive - University of Gothenburg/ARMS MBON/data paper/ARMS_data_paper_Nauras/Analysis/Nauras/COI")

# read Extended_final_tables of each sequencing run
# One Gdynia sample was processed separately and therefore has its own table
# separate ASV_xy prefixes from long ID code in first column

ASVcounts_April2021 <- read_xlsx("Extended_final_table_April2021_COI_noBlank.xlsx")
ASVcounts_April2021 <- separate(data = ASVcounts_April2021, col = "ASV_number:amplicon", into = c("ASV_number", "ID"), sep = "\\:")

ASVcounts_Jan2020 <- read_xlsx("Extended_final_table_January2020_COI_noBlank.xlsx")
ASVcounts_Jan2020 <- separate(data = ASVcounts_Jan2020, col = "ASV_number:amplicon", into = c("ASV_number", "ID"), sep = "\\:")

ASVcounts_Jan2022 <- read_xlsx("Extended_final_table_January2022_COI_noBlank.xlsx")
ASVcounts_Jan2022 <- separate(data = ASVcounts_Jan2022, col = "ASV_number:amplicon", into = c("ASV_number", "ID"), sep = "\\:")

ASVcounts_July2019 <- read_xlsx("Extended_final_table_July2019_COI_noBlank.xlsx")
ASVcounts_July2019 <- separate(data = ASVcounts_July2019, col = "ASV_number:amplicon", into = c("ASV_number", "ID"), sep = "\\:")

ASVcounts_May2021 <- read_xlsx("Extended_final_table_May2021_COI_noBlank.xlsx")
ASVcounts_May2021 <- separate(data = ASVcounts_May2021, col = "ASV_number:amplicon", into = c("ASV_number", "ID"), sep = "\\:")

ASVcounts_Sep2020 <- read_xlsx("Extended_final_table_September2020_COI_noBlank.xlsx")
ASVcounts_Sep2020 <- separate(data = ASVcounts_Sep2020, col = "ASV_number:amplicon", into = c("ASV_number", "ID"), sep = "\\:")

ASVcounts_Aug2023 <- read_xlsx("Extended_final_table_August2023_COI_noBlank.xlsx")
ASVcounts_Aug2023 <- separate(data = ASVcounts_Aug2023, col = "ASV_number:amplicon", into = c("ASV_number", "ID"), sep = "\\:")
colnames(ASVcounts_Aug2023)[3:(ncol(ASVcounts_Aug2023)-2)]<-paste0(colnames(ASVcounts_Aug2023)[3:(ncol(ASVcounts_Aug2023)-2)], "_rep") # Because the samples of August2023 are re-sequenced samples, we need to add a string to them to distinguish their names from the samples with the same name

ASVcounts_GDY <- read_xlsx("Extended_final_table_ARMS_Gdynia_GDY1_20180813_20191029_SF38_DMSO_COI_noBlank.xlsx")
ASVcounts_GDY <- separate(data = ASVcounts_GDY, col = "ASV_number:amplicon", into = c("ASV_number", "ID"), sep = "\\:")

# read tax_assignment files of each sequencing run 

Taxfile_April2021 <- read.table("tax_assignments_April2021_COI_noBlank.tsv",sep="\t")
Taxfile_Jan2020 <- read.table("tax_assignments_January2020_COI_noBlank.tsv",sep="\t")
Taxfile_Jan2022 <- read.table("tax_assignments_January2022_COI_noBlank.tsv",sep="\t")
Taxfile_July2019 <- read.table("tax_assignments_July2019_COI_noBlank.tsv",sep="\t")
Taxfile_May2021 <- read.table("tax_assignments_May2021_COI_noBlank.tsv",sep="\t")
Taxfile_Sep2020 <- read.table("tax_assignments_September2020_COI_noBlank.tsv",sep="\t")
Taxfile_Aug2023 <- read.table("tax_assignments_August2023_COI_noBlank.tsv",sep="\t")
Taxfile_GDY <- read.table("tax_assignments_ARMS_Gdynia_GDY1_20180813_20191029_SF38_DMSO_COI_noBlank.tsv",sep="\t")

# Make lists of Taxfiles and ASVcounts

all_Taxfile<-tibble::lst(Taxfile_April2021,Taxfile_Jan2020,Taxfile_Jan2022,Taxfile_July2019,Taxfile_May2021,Taxfile_Sep2020,Taxfile_Aug2023,Taxfile_GDY) # lst from tibble will preserve object names
all_ASVcounts<-tibble::lst(ASVcounts_April2021,ASVcounts_Jan2020,ASVcounts_Jan2022,ASVcounts_July2019,ASVcounts_May2021,ASVcounts_Sep2020,ASVcounts_Aug2023,ASVcounts_GDY) # lst from tibble will preserve object names

# set taxonomy assignments with confidence below 0.8 to NA and subset to ASVs in ASVcounts

for(i in 1:length(all_Taxfile)) { 
  colnames(all_Taxfile[[i]])[1]<-"ID"
  all_Taxfile[[i]] <- separate(data = all_Taxfile[[i]], col = ID, into = c("ID", "Read_number"), sep = "\\_")
  all_Taxfile[[i]] <-all_Taxfile[[i]][,-c(2,4,7,10,13,16,19,22)]
  all_Taxfile[[i]][,2]<-ifelse(all_Taxfile[[i]][,3]<=0.8,NA,all_Taxfile[[i]][,2])
  all_Taxfile[[i]][,4]<-ifelse(all_Taxfile[[i]][,5]<=0.8,NA,all_Taxfile[[i]][,4])
  all_Taxfile[[i]][,6]<-ifelse(all_Taxfile[[i]][,7]<=0.8,NA,all_Taxfile[[i]][,6])
  all_Taxfile[[i]][,8]<-ifelse(all_Taxfile[[i]][,9]<=0.8,NA,all_Taxfile[[i]][,8])
  all_Taxfile[[i]][,10]<-ifelse(all_Taxfile[[i]][,11]<=0.8,NA,all_Taxfile[[i]][,10])
  all_Taxfile[[i]][,12]<-ifelse(all_Taxfile[[i]][,13]<=0.8,NA,all_Taxfile[[i]][,12])
  all_Taxfile[[i]][,14]<-ifelse(all_Taxfile[[i]][,15]<=0.8,NA,all_Taxfile[[i]][,14])
  all_Taxfile[[i]]<-all_Taxfile[[i]][,-c(3,5,7,9,11,13,15)]
  colnames(all_Taxfile[[i]])<-c("ID","Kingdom","Phylum","Class","Order","Family","Genus","Species")
  all_Taxfile[[i]] <- all_Taxfile[[i]] %>% filter(ID %in% all_ASVcounts[[i]]$ID) # Just precautionary measure in case taxonomy files still contain ASVs that were removed during blank-correction
}  

# read fasta files to get sequences 
# all_sequences_grouped_XX_COI.fa files contain all clustered ASV sequences, including ASVs present in blank samples and ASVs not classified taxonomically
# the following fasta files therefore contain more sequences than ASVs present in the Extended_final_table_XX_COI_noBlank.xlsx files

fasta_April2021 <- readDNAStringSet("all_sequences_grouped_April2021_COI.fa")

fasta_Jan2020 <- readDNAStringSet("all_sequences_grouped_January2020_COI.fa")

fasta_Jan2022 <- readDNAStringSet("all_sequences_grouped_January2022_COI.fa")

fasta_July2019 <- readDNAStringSet("all_sequences_grouped_July2019_COI.fa")

fasta_May2021 <- readDNAStringSet("all_sequences_grouped_May2021_COI.fa")

fasta_Sep2020 <- readDNAStringSet("all_sequences_grouped_September2020_COI.fa")

fasta_Aug2023 <- readDNAStringSet("all_sequences_grouped_August2023_COI.fa")

fasta_GDY <- readDNAStringSet("all_sequences_grouped_ARMS_Gdynia_GDY1_20180813_20191029_SF38_DMSO_COI.fa")

# Make list of fastas

all_fasta<-tibble::lst(fasta_April2021,fasta_Jan2020,fasta_Jan2022,fasta_July2019,fasta_May2021,fasta_Sep2020,fasta_Aug2023,fasta_GDY) 

# Make fasta files into dataframes and remove read counts from sequence identifiers

IDs<-list()
sequences<-list()
all_seqs<-list()
for(i in 1:length(all_fasta)) {
  IDs[[i]]<-names(all_fasta[[i]])
  sequences[[i]]<-paste(all_fasta[[i]])
  all_seqs[[i]]<-data.frame(IDs[[i]],sequences[[i]])
  colnames(all_seqs[[i]])<-c("ID","seqs")
  all_seqs[[i]]$ID<-gsub("\\_.*","",all_seqs[[i]]$ID)
}  

# Get number of unique ASVs prior to any curation

length(unique(rbind(all_seqs[[1]],all_seqs[[2]],all_seqs[[3]],all_seqs[[4]],all_seqs[[5]],all_seqs[[6]],all_seqs[[7]],all_seqs[[8]])$seqs))

# Subset to ASVs found in count tables (these are the ASVs with a taxonomic classification remaining after blank-curation)

for(i in 1:length(all_seqs)) {
  all_seqs[[i]] <- all_seqs[[i]] %>% filter(ID %in% all_ASVcounts[[i]]$ID)
}  

# Add sequences to Taxfiles

for(i in 1:length(all_Taxfile)) {
  for(j in 1:nrow(all_Taxfile[[i]])) {
    all_Taxfile[[i]]$seqs[all_Taxfile[[i]]$ID == all_seqs[[i]]$ID[j]] <- all_seqs[[i]][,2][j]
  }
}

# Add sequences to ASVcounts

for(i in 1:length(all_ASVcounts)) {
  for(j in 1:nrow(all_ASVcounts[[i]])) {
    all_ASVcounts[[i]]<-data.frame(all_ASVcounts[[i]])
    all_ASVcounts[[i]]$seqs[all_ASVcounts[[i]]$ID == all_seqs[[i]]$ID[j]] <- all_seqs[[i]][,2][j]
  }
}

# Bring sequence columns to front, remove ID and other columns

all_Taxfile<-lapply(all_Taxfile,function(x) x %>% relocate(seqs))
all_Taxfile<-lapply(all_Taxfile,function(x) x %>% select(-ID))

all_ASVcounts<-lapply(all_ASVcounts,function(x) x %>% relocate(seqs))
all_ASVcounts<-lapply(all_ASVcounts,function(x) x %>% select(-c(ASV_number,ID,Classification,TAXON.NCBI_TAX_ID)))

# Combine Taxfiles and ASVcounts among each other, set NAs to zero

ASVcounts<-rbind.fill(all_ASVcounts[[1]],all_ASVcounts[[2]],all_ASVcounts[[3]],all_ASVcounts[[4]],all_ASVcounts[[5]],all_ASVcounts[[6]],all_ASVcounts[[7]],all_ASVcounts[[8]])
ASVcounts[is.na(ASVcounts)]<-0

Taxfile<-rbind.fill(all_Taxfile[[1]],all_Taxfile[[2]],all_Taxfile[[3]],all_Taxfile[[4]],all_Taxfile[[5]],all_Taxfile[[6]],all_Taxfile[[7]],all_Taxfile[[8]])

# Aggregate ASV counts based on unique sequences

ASVcounts<-aggregate(.~ seqs,data = ASVcounts,FUN=sum)

# Keep only distinct sequence-taxonomy combinations

Taxfile<-Taxfile %>% distinct()

# Sort Taxfile based on order in ASVcounts

Taxfile<-Taxfile[order(match(Taxfile[,1],ASVcounts[,1])),]

# In the Species column

Taxfile$Species <- ifelse(grepl("*_sp$", Taxfile$Species), NA, Taxfile$Species) # sets all Species entries that end with just _sp as NA

# Make fasta file of unique sequences with short ASV names and save to file

fasta<-data.frame(ASV=paste0(">ASV",1:nrow(ASVcounts)),seqs=ASVcounts$seqs) |>
  pivot_longer(everything()) |> 
  subset(select=-name)

write.table(fasta,"COI_ASVs.fa",row.names = F,col.names = F,quote=F)

# Write ASV table to file

ASVcounts$ASV<-paste0("ASV",1:nrow(ASVcounts))
ASVcounts<-ASVcounts %>% relocate(ASV) %>% select(-seqs)
write.table(ASVcounts,"COI_ASV_table.txt",sep="\t",row.names = F)

# Write taxonomy table to file

Taxfile$ASV<-paste0("ASV",1:nrow(Taxfile))
Taxfile<-Taxfile %>% relocate(ASV) %>% select(-seqs)
write.table(Taxfile,"COI_tax_table.txt",sep="\t",row.names = F)

# get number of unique, classified, blank-corrected ASVs

length(ASVcounts$ASV)

# Get number of reads

sum(colSums(ASVcounts[,-1]))
