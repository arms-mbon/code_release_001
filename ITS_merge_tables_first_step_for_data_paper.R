library(tidyr)
library(tibble)
library(Biostrings)
library(dplyr)                                                 
library(plyr)                                                  
library(readxl)
library(stringr)

setwd("C:/Users/Nauras/OneDrive - University of Gothenburg/ARMS MBON/data paper/ARMS_data_paper_Nauras/Analysis/Nauras/ITS")

# read Extended_final_tables of each sequencing run

ASVcounts_April2021 <- read_xlsx("Extended_final_table_April2021_ITS_noBlank.xlsx")

ASVcounts_July2019 <- read_xlsx("Extended_final_table_July2019_ITS_noBlank.xlsx")

# read fasta files to get sequences
# all_sequences_grouped_XX_ITS.fa files contain all clustered ASV sequences, including ASVs present in blank samples andASVs not classified taxonomically
# the following fasta files therefore contain more sequences thanASVs present in the Extended_final_table_XX_ITS_noBlank.xlsx files

fasta_April2021 <- readDNAStringSet("all_sequences_grouped_April2021_ITS.fa")

fasta_July2019 <- readDNAStringSet("all_sequences_grouped_July2019_ITS.fa")

# Make lists of ASVcounts and fastas

all_ASVcounts<-tibble::lst(ASVcounts_April2021,ASVcounts_July2019) # lst from tibble will preserve object names
all_ASVcounts<-lapply(all_ASVcounts,function(x) {colnames(x)[1]<-"ASV";x})
all_fasta<-tibble::lst(fasta_April2021,fasta_July2019) 

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

# Get number of unique ASVs prior to any curation

length(unique(rbind(all_seqs[[1]],all_seqs[[2]])$seqs))

# Subset to ASVs found in count tables (these are the ASVs with a taxonomic classification remaining after blank-curation)

for(i in 1:length(all_seqs)) {
  all_seqs[[i]] <- all_seqs[[i]] %>% filter(ID %in% all_ASVcounts[[i]]$ASV)
}  

# Add sequences to ASVcounts

for(i in 1:length(all_ASVcounts)) {
  for(j in 1:nrow(all_ASVcounts[[i]])) {
    all_ASVcounts[[i]]<-data.frame(all_ASVcounts[[i]])
    all_ASVcounts[[i]]$seqs[all_ASVcounts[[i]]$ASV == all_seqs[[i]]$ID[j]] <- all_seqs[[i]][,2][j]
  }
}

# Make taxonomy tables

taxa<-lapply(all_ASVcounts,function(x) x %>% select(c(seqs,Classification)))

# Bring sequence columns to front, remove ASV and taxonomy columns

all_ASVcounts<-lapply(all_ASVcounts,function(x) x %>% relocate(seqs))
all_ASVcounts<-lapply(all_ASVcounts,function(x) x %>% select(-c(ASV,Classification,TAXON.NCBI_TAX_ID)))

# Combine ASVcounts, set NAs to zero

ASVcounts<-rbind.fill(all_ASVcounts[[1]],all_ASVcounts[[2]])
ASVcounts[is.na(ASVcounts)]<-0

# Aggregate ASV counts based on unique sequences

ASVcounts<-aggregate(.~ seqs,data = ASVcounts,FUN=sum)

# Combine taxonomy tables

Taxfile<-rbind(taxa[[1]],taxa[[2]])

# Bring taxonomy strings into separate columns and set certain entries as NA

Taxfile<-separate_wider_delim(Taxfile,Classification, delim = ";", names_sep="",too_few = "align_start")

Taxfile[] <- lapply(Taxfile, function(x) replace(x, grepl(" |_sp$|Incertae", x), NA)) # sets all levels including a space or Incertae Sedis or levels ending in _sp as NA

Taxfile<-Taxfile[colSums(!is.na(Taxfile)) > 0] # Remove columns which are now left with NAs only

# Keep only distinct sequence-taxonomy combinations

Taxfile<-Taxfile %>% distinct()

# Sort Taxfile based on order in ASVcounts

Taxfile<-as.data.frame(Taxfile)
Taxfile<-Taxfile[order(match(Taxfile[,1],ASVcounts[,1])),]

# Set taxonomic levels

colnames(Taxfile)[2:ncol(Taxfile)]<-c("Domain","Supergroup","Kingdom","Phylum","Class","Order","Family","Genus","Species")

# Make fasta file of unique sequences with short ASV names and save to file

fasta<-data.frame(ASV=paste0(">ASV",1:nrow(ASVcounts)),seqs=ASVcounts$seqs) |>
  pivot_longer(everything()) |> 
  subset(select=-name)

write.table(fasta,"ITS_ASVs.fa",row.names = F,col.names = F,quote=F)

# Write ASV table to file

ASVcounts$ASV<-paste0("ASV",1:nrow(ASVcounts))
ASVcounts<-ASVcounts %>% relocate(ASV) %>% select(-seqs)
write.table(ASVcounts,"ITS_ASV_table.txt",sep="\t",row.names = F)

# Write taxonomy table to file

Taxfile$ASV<-paste0("ASV",1:nrow(Taxfile))
Taxfile<-Taxfile %>% relocate(ASV) %>% select(-seqs)
write.table(Taxfile,"ITS_tax_table.txt",sep="\t",row.names = F)

# get number of unique, classified, blank-corrected ASVs

length(ASVcounts$ASV)

# Get number of reads

sum(colSums(ASVcounts[,-1]))
