library(dplyr)
library(ggplot2)
library(egg)
library(grafify)
library(phyloseq)
library(vegan)
library(grid) # had to be called, even though it should come with ggplot??
library(ggpubr)
library(scales)
library(tidyr)
library(xlsx)
library(readxl)
library(writexl)
library(UpSetR)

setwd("~/ARMS_data_paper_Nauras/Analysis")

### Create phyloseq objects, do first quick data assessments and get info on recovered phyla species-level assignments ###

## COI ##

# Read ASV counts

ASVcountsCOI<-read.table("COI/COI_ASV_table.txt",header=T,check.names=F, sep="\t",row.names = 1)

# Read ASV taxonomy (needs to be read as matrix, may cause problems otherwise when phyloseq object will be created)

ASVtaxaCOI<-as.matrix(read.table("COI/COI_tax_table.txt",header=T,check.names=F, sep="\t",row.names = 1))

# Sort count table based on order in tax table (precautionary measure)

ASVcountsCOI<-ASVcountsCOI[order(match(rownames(ASVcountsCOI),rownames(ASVtaxaCOI))),]

# Read sample metadata and format a bit

samples<-read.table("sample_data.txt",header=T,check.names=F,sep="\t") 
rownames(samples) <- samples$Sample
colnames(samples)[1]<-"MaterialSampleID"
samples$MaterialSampleID<-gsub("_r1","",samples$MaterialSampleID) # Remove the "_r1" strings in MaterialSampleID column
samples$MaterialSampleID<-gsub("_r2","",samples$MaterialSampleID) # Remove the "_r2" strings in MaterialSampleID column

# Create phyloseq object

psCOI <- phyloseq(otu_table(ASVcountsCOI,taxa_are_rows = TRUE), sample_data(samples), tax_table(ASVtaxaCOI))

# Get number of samples that produced ASVs through PEMA processing

psCOI

# Remove the sediment and plankton samples (some sediment and plankton samples were sequenced as a trial during the initial phase of the ARMS program)

psCOI <- subset_samples(psCOI,Fraction!="SED")
psCOI <- subset_samples(psCOI,Fraction!="PS")

# Remove certain erroneous ASVs

to_remove_taxa_COI<-c("Homo","Bos","Canis","Gynaikothrips","Dorypteryx","Fannia","Bactrocera")

erroneous_coi<-subset_taxa(psCOI,Genus %in% to_remove_taxa_COI) # phyloseq object with erroneous ASVs
erroneous_coi <- prune_samples(sample_sums(erroneous_coi) > 0, erroneous_coi) # Remove samples with are left with a read number of zero
xlsx::write.xlsx(as.data.frame(cbind(otu_table(erroneous_coi),tax_table(erroneous_coi))), file = "erroneus_sequences_removed.xlsx",sheetName = "COI", append = FALSE) # Save to file

psCOI<-subset_taxa(psCOI,!Genus %in% to_remove_taxa_COI) # Make phyloseq object without these sequences

# Remove samples with are left with a read number of zero

psCOI <- prune_samples(sample_sums(psCOI) > 0, psCOI)

# Some samples were re-sequenced in August 2023 and in most of those cases two sequence samples are therefore present in the data set
# Assess rarefaction curves, ASV counts and taxonomy for each of those sample pairs and remove sample with lower diversity or taxonomic resolution

katza<-subset_samples(psCOI,MaterialSampleID=="ARMS_Eilat_Katza1_20181024_20200706_MF500")
katza<-prune_taxa(rowSums(otu_table(katza))>0,katza)
rarecurve(t(otu_table(katza)), step=50, cex=0.5)
katza<-cbind(tax_table(katza),otu_table(katza))
katza

torallaA<-subset_samples(psCOI,MaterialSampleID=="ARMS_Vigo_TorallaA_20190625_20191014_MF500")
torallaA<-prune_taxa(rowSums(otu_table(torallaA))>0,torallaA)
rarecurve(t(otu_table(torallaA)), step=50, cex=0.5)
torallaA<-cbind(tax_table(torallaA),otu_table(torallaA))
torallaA

torallaB<-subset_samples(psCOI,MaterialSampleID=="ARMS_Vigo_TorallaB_20190625_20191014_MF500")
torallaB<-prune_taxa(rowSums(otu_table(torallaB))>0,torallaB)
rarecurve(t(otu_table(torallaB)), step=50, cex=0.5)
torallaB<-cbind(tax_table(torallaB),otu_table(torallaB))
torallaB

torallaC<-subset_samples(psCOI,MaterialSampleID=="ARMS_Vigo_TorallaC_20190625_20191014_MF500")
torallaC<-prune_taxa(rowSums(otu_table(torallaC))>0,torallaC)
rarecurve(t(otu_table(torallaC)), step=50, cex=0.5)
torallaC<-cbind(tax_table(torallaC),otu_table(torallaC))
torallaC

gdynia<-subset_samples(psCOI,MaterialSampleID=="ARMS_Gdynia_GDY1_20180813_20191029_SF38_DMSO")
gdynia<-prune_taxa(rowSums(otu_table(gdynia))>0,gdynia)
rarecurve(t(otu_table(gdynia)), step=50, cex=0.5)
gdynia<-cbind(tax_table(gdynia),otu_table(gdynia))
gdynia

fornace<-subset_samples(psCOI,MaterialSampleID=="ARMS_GulfOfPiran_Fornace_20180815_20181118_SF_ETOH")
fornace<-prune_taxa(rowSums(otu_table(fornace))>0,fornace)
rarecurve(t(otu_table(fornace)), step=50, cex=0.5)
fornace<-cbind(tax_table(fornace),otu_table(fornace))
fornace

to_remove_samples<-c("ARMS_Eilat_Katza1_20181024_20200706_MF500_r1",
                     "ARMS_Vigo_TorallaA_20190625_20191014_MF500_r1",
                     "ARMS_Vigo_TorallaB_20190625_20191014_MF500_r1",
                     "ARMS_Vigo_TorallaC_20190625_20191014_MF500_r1",
                     "ARMS_Gdynia_GDY1_20180813_20191029_SF38_DMSO_r1",
                     "ARMS_GulfOfPiran_Fornace_20180815_20181118_SF_ETOH_r1")

psCOI <- prune_samples(!(sample_names(psCOI) %in% to_remove_samples), psCOI)

# Some samples have been preserved in DMSO as well as EtOH initially as a trial. 
# Where for the remaining samples both replicates exist, keep only DMSO samples.

sort(sample_names(psCOI)) # Check sample names to see where both replicates exist (check for same sample name with and without "_EtOH"

to_remove_pres_COI <-c("ARMS_Crete_1HERP_20180928_20190128_MF125_ETOH_r1",
                   "ARMS_Crete_1HERP_20180928_20190128_MF500_ETOH_r1",
                   "ARMS_Crete_1HERP_20180928_20190128_SF_ETOH_r1",
                   "ARMS_GulfOfPiran_Fornace_20180815_20181118_MF100_ETOH_r1",
                   "ARMS_GulfOfPiran_Fornace_20180815_20181118_MF500_ETOH_r2",
                   "ARMS_GulfOfPiran_Fornace_20180815_20181118_SF_ETOH_r2",
                   "ARMS_Koster_VH2_20180418_20180906_SF40_ETOH_r1",
                   "ARMS_Plymouth_MBA1_20180701_20181001_MF100_ETOH_r1",
                   "ARMS_Plymouth_MBA1_20180701_20181001_MF500_ETOH_r1",
                   "ARMS_Plymouth_MBA1_20180701_20181001_SF40_ETOH_r1",
                   "ARMS_Plymouth_MBA2_20180701_20181001_SF40_ETOH_r1",
                   "ARMS_Roscoff_BasBloS1_20180711_20181025_SF40_ETOH_A_r1",
                   "ARMS_Roscoff_BasBloS1_20180711_20181025_SF40_ETOH_B_r1",
                   "ARMS_Roscoff_MarBloR1_20180709_20181024_MF100_ETOH_r1",
                   "ARMS_Roscoff_MarBloR1_20180709_20181024_MF500_ETOH_r1",
                   "ARMS_Roscoff_MarBloR1_20180709_20181024_SF40_ETOH_r1")

psCOI <- prune_samples(!(sample_names(psCOI) %in% to_remove_pres_COI), psCOI)

# Two samples from Roscoff, France have been processed as replicates. #

# Assess rarefaction curves to see which sample to keep
ps_roscoff_COI<-prune_samples(sample_names(psCOI) %in% c("ARMS_Roscoff_BasBloS1_20180711_20181025_SF40_DMSO_A_r1","ARMS_Roscoff_BasBloS1_20180711_20181025_SF40_DMSO_B_r1"), psCOI)
rarecurve(t(otu_table(ps_roscoff_COI)), step=50, cex=0.5)

# Remove the sample with lower sample size/ASV richness
psCOI <- prune_samples(!(sample_names(psCOI) %in% "ARMS_Roscoff_BasBloS1_20180711_20181025_SF40_DMSO_B_r1"), psCOI)

# Remove ASVs which have a total abundance of zero after removing samples during the previous step

psCOI<-prune_taxa(rowSums(otu_table(psCOI))>0,psCOI)

# Check number of ASVs after further curation and filtering

psCOI

# Get total number of reads

sum(sample_sums(psCOI))

# Get number of ASVs assigned to species level

subset_taxa(psCOI,!is.na(Species))

# Get number of unique species identified with Linnean name 

length(unique(tax_table(subset_taxa(psCOI,!is.na(Species)))[,ncol(tax_table(psCOI))]))

# Get number of species observations for occurrences with a minimum of 2 reads (i.e. all presence-absence occurrences of all ASVs classified to species level across all samples)

psCOI_species<-subset_taxa(psCOI,!is.na(Species))

psCOI_species_obs<-psCOI_species

otu_table(psCOI_species_obs)[otu_table(psCOI_species_obs)<2]<-0

otu_table(psCOI_species_obs)[otu_table(psCOI_species_obs)>1]<-1

sum(sample_sums(psCOI_species_obs))

# Get data set agglomerated at phylum level and with relative abundances #

psCOI.phylum <- tax_glom(psCOI, taxrank = "Phylum",NArm=F) 

phylum_COI<-as.data.frame(cbind(as.data.frame(tax_table(psCOI.phylum))[,2],taxa_sums(psCOI.phylum)))

colnames(phylum_COI)<-c("Phylum","reads")

phylum_COI$reads<-as.numeric(phylum_COI$reads)

phylum_COI[is.na(phylum_COI)]<-"NA"

phylum_COI<-aggregate(.~ Phylum,data = phylum_COI,FUN=sum)

phylum_COI<-phylum_COI[order(phylum_COI$reads, decreasing = TRUE),]  

phylum_COI$reads<-phylum_COI$reads/sum(phylum_COI$reads)

# Change classification in the phylum level column to correct phylum name
# This was done based on web-based search using WoRMS or other scientific publications

phylum_COI$Phylum # check all unique entries in the Phylum column of the COI data set

phylum_COI$Phylum<-ifelse(grepl("Ectocarpales|Fucales|Eustigmatales|Laminariales|Cutleriales|Dictyotales|Chattonellales|Chromulinales|Parmales|Sphacelariales|Florenciellales|Goniochloridales", phylum_COI$Phylum),"Ochrophyta",phylum_COI$Phylum)
phylum_COI$Phylum<-ifelse(grepl("Discosea|Evosea", phylum_COI$Phylum),"Amoebozoa",phylum_COI$Phylum)
phylum_COI$Phylum<-ifelse(grepl("Peridiniales|Syndiniales|Suessiales|Gymnodiniales|Gonyaulacales", phylum_COI$Phylum),"Myzozoa",phylum_COI$Phylum)
phylum_COI$Phylum<-ifelse(grepl("Myzocytiopsidales|Saprolegniales|Anisolpidiales|Peronosporales|Pythiales", phylum_COI$Phylum),"Oomycota",phylum_COI$Phylum)
phylum_COI$Phylum<-ifelse(grepl("Imbricatea", phylum_COI$Phylum),"Cercozoa",phylum_COI$Phylum)
phylum_COI$Phylum<-ifelse(grepl("Apusomonadidae", phylum_COI$Phylum),"Apusozoa",phylum_COI$Phylum)
phylum_COI$Phylum<-ifelse(grepl("Ministeria", phylum_COI$Phylum),"Choanozoa",phylum_COI$Phylum)
phylum_COI$Phylum<-ifelse(grepl("Thraustochytrida|Bicosoecida", phylum_COI$Phylum),"Bigyra",phylum_COI$Phylum)

phylum_COI<-aggregate(.~ Phylum,data = phylum_COI,FUN=sum)

phylum_COI<-phylum_COI[order(phylum_COI$reads, decreasing = TRUE),]  

# get number and names of identified phyla
phyla_identified_COI<-phylum_COI
length(which(phyla_identified_COI$Phylum!="NA")) # Number of recovered phyla
phyla_identified_COI[phyla_identified_COI$Phylum=="NA",] # Percentage of reads unclassified at phylum level
write.table(phyla_identified_COI,"phyla_identified_COI.txt",sep="\t",row.names = F,quote = F)

phylum_COI$Phylum[12:nrow(phylum_COI)]<-"Other" # keep only top 10 most abundant phyla and NA and set rest to "Other"

phylum_COI<-aggregate(.~ Phylum,data = phylum_COI,FUN=sum)

phylum_COI<-phylum_COI[order(phylum_COI$reads, decreasing = TRUE),]

phylum_COI<-phylum_COI[c(2:10,12,11,1),] # Re-order entries for plot later on

# For ASVs with species level classification, get data set agglomerated at phylum  level 

taxa_COI<-as.data.frame(tax_table(psCOI_species))

specs_COI<-distinct(taxa_COI,Species,.keep_all=T)

specs_COI$abundance<-rep(1,nrow(specs_COI))

specs_COI<-specs_COI %>% select(c(Phylum,abundance))

specs_COI<-aggregate(.~ Phylum,data = specs_COI,FUN=sum)

specs_COI<-specs_COI[order(specs_COI$abundance, decreasing = TRUE),]

# Change classification in the phylum level column to correct phylum name
# This was done based on web-based search using WoRMS or other scientific publications
# We just use the code from above here, because the species-assigned data set is just a subset of the data set above

specs_COI$Phylum<-ifelse(grepl("Ectocarpales|Fucales|Eustigmatales|Laminariales|Cutleriales|Dictyotales|Chattonellales|Chromulinales|Parmales|Sphacelariales|Florenciellales|Goniochloridales", specs_COI$Phylum),"Ochrophyta",specs_COI$Phylum)
specs_COI$Phylum<-ifelse(grepl("Discosea|Evosea", specs_COI$Phylum),"Amoebozoa",specs_COI$Phylum)
specs_COI$Phylum<-ifelse(grepl("Peridiniales|Syndiniales|Suessiales|Gymnodiniales|Gonyaulacales", specs_COI$Phylum),"Myzozoa",specs_COI$Phylum)
specs_COI$Phylum<-ifelse(grepl("Myzocytiopsidales|Saprolegniales|Anisolpidiales|Peronosporales|Pythiales", specs_COI$Phylum),"Oomycota",specs_COI$Phylum)
specs_COI$Phylum<-ifelse(grepl("Imbricatea", specs_COI$Phylum),"Cercozoa",specs_COI$Phylum)
specs_COI$Phylum<-ifelse(grepl("Apusomonadidae", specs_COI$Phylum),"Apusozoa",specs_COI$Phylum)
specs_COI$Phylum<-ifelse(grepl("Ministeria", specs_COI$Phylum),"Choanozoa",specs_COI$Phylum)
specs_COI$Phylum<-ifelse(grepl("Thraustochytrida|Bicosoecida", specs_COI$Phylum),"Bigyra",specs_COI$Phylum)

specs_COI<-aggregate(.~ Phylum,data = specs_COI,FUN=sum)

specs_COI<-specs_COI[order(specs_COI$abundance, decreasing = TRUE),]  

# get number and names of identified phyla
phyla_specs_identified_COI<-specs_COI
length(phyla_specs_identified_COI$Phylum) # Number of recovered phyla
write.table(phyla_specs_identified_COI,"phyla_specs_identified_COI.txt",sep="\t",row.names = F,quote = F)

specs_COI$Phylum<-ifelse(specs_COI$abundance>=10,specs_COI$Phylum,"Other") # keep only phyla with at least 10 species and set rest to "Other"

specs_COI<-aggregate(.~ Phylum,data = specs_COI,FUN=sum)

specs_COI<-specs_COI[order(specs_COI$abundance, decreasing = TRUE),]

specs_COI<-specs_COI[c(1:5,7:15,6),] # re-order for plot later on

## 18S ##

# Read OTU counts

OTUcounts18S<-read.table("18S/18S_OTU_table.txt",header=T,check.names=F, sep="\t",row.names = 1)

# Read OTU taxonomy (needs to be read as matrix, may cause problems otherwise when phyloseq object will be created)

OTUtaxa18S<-as.matrix(read.table("18S/18S_tax_table.txt",header=T,check.names=F, sep="\t",row.names = 1))

# Sort count table based on order in tax table (precautionary measure)

OTUcounts18S<-OTUcounts18S[order(match(rownames(OTUcounts18S),rownames(OTUtaxa18S))),]

# Create phyloseq object

ps18S <- phyloseq(otu_table(OTUcounts18S,taxa_are_rows = TRUE), sample_data(samples), tax_table(OTUtaxa18S))

# Get number of samples that produced OTUs through PEMA processing

ps18S

# Remove the sediment and plankton samples (some sediment and plankton samples were sequenced as a trial during the initial phase of the ARMS program)

ps18S <- subset_samples(ps18S,Fraction!="SED")
ps18S <- subset_samples(ps18S,Fraction!="PS")

# Remove certain erroneous OTUs

to_remove_taxa_18S<-c("Zea_mays","Stegobium_paniceum")

erroneous_18S_1<-subset_taxa(ps18S,Species %in% to_remove_taxa_18S) # phyloseq object 1 with erroneous OTUs
erroneous_18S_2<-subset_taxa(ps18S,Level_XXX %in% "Drosophila") # phyloseq object 2 with erroneous OTUs
erroneous_18S<-merge_phyloseq(erroneous_18S_1,erroneous_18S_2) # merge both objects
erroneous_18S <- prune_samples(sample_sums(erroneous_18S) > 0, erroneous_18S) # Remove samples with are left with a read number of zero
xlsx::write.xlsx(as.data.frame(cbind(otu_table(erroneous_18S),tax_table(erroneous_18S))), file = "erroneus_sequences_removed.xlsx",sheetName = "18S", append = T) # Save to file

ps18S<-subset_taxa(ps18S,!Species %in% to_remove_taxa_18S) # Make phyloseq object without these sequences
ps18S<-subset_taxa(ps18S,!Level_XXX %in% "Drosophila") # Make phyloseq object without these sequences

# Remove samples with are left with a read number of zero

ps18S <- prune_samples(sample_sums(ps18S) > 0, ps18S)

# Some samples were re-sequenced in August 2023 and in most of those cases two sequence samples are therefore present in the data set
# Assess rarefaction curves, OTU counts and taxonomy for each of those sample pairs and remove sample with lower diversity or taxonomic resolution

gdynia<-subset_samples(ps18S,MaterialSampleID=="ARMS_Gdynia_GDY1_20180813_20191029_SF38_DMSO")
gdynia<-prune_taxa(rowSums(otu_table(gdynia))>0,gdynia)
rarecurve(t(otu_table(gdynia)), step=50, cex=0.5)
gdynia<-cbind(tax_table(gdynia),otu_table(gdynia))
gdynia

koster<-subset_samples(ps18S,MaterialSampleID=="ARMS_Koster_VH1_20190527_20200716_MF100")
koster<-prune_taxa(rowSums(otu_table(koster))>0,koster)
rarecurve(t(otu_table(koster)), step=50, cex=0.5)
koster<-cbind(tax_table(koster),otu_table(koster))
koster

fornace<-subset_samples(ps18S,MaterialSampleID=="ARMS_GulfOfPiran_Fornace_20180815_20181118_MF500_ETOH")
fornace<-prune_taxa(rowSums(otu_table(fornace))>0,fornace)
rarecurve(t(otu_table(fornace)), step=50, cex=0.5)
fornace<-cbind(tax_table(fornace),otu_table(fornace))
fornace

to_remove_samples<-c("ARMS_Gdynia_GDY1_20180813_20191029_SF38_DMSO_r1","ARMS_Koster_VH1_20190527_20200716_MF100_r2","ARMS_GulfOfPiran_Fornace_20180815_20181118_MF500_ETOH_r2")

ps18S <- prune_samples(!(sample_names(ps18S) %in% to_remove_samples), ps18S)

# Some samples have been preserved in DMSO as well as EtOH initially as a trial. 
# Where for the remaining samples both replicates exist, keep only DMSO samples.

sort(sample_names(ps18S)) # Check sample names to see where both replicates exist (check for same sample name with and without "_EtOH"

to_remove_pres_18S <-c("ARMS_Crete_1HERP_20180928_20190128_MF125_ETOH_r1",
                       "ARMS_Crete_1HERP_20180928_20190128_MF500_ETOH_r1",
                       "ARMS_Crete_1HERP_20180928_20190128_SF_ETOH_r1",
                       "ARMS_GulfOfPiran_Fornace_20180815_20181118_MF100_ETOH_r1",
                       "ARMS_GulfOfPiran_Fornace_20180815_20181118_MF500_ETOH_r1",
                       "ARMS_GulfOfPiran_Fornace_20180815_20181118_SF_ETOH_r2",
                       "ARMS_Koster_VH2_20180418_20180906_SF40_ETOH_r1",
                       "ARMS_Plymouth_MBA1_20180701_20181001_MF100_ETOH_r1",
                       "ARMS_Plymouth_MBA1_20180701_20181001_MF500_ETOH_r1",
                       "ARMS_Plymouth_MBA1_20180701_20181001_SF40_ETOH_r1",
                       "ARMS_Plymouth_MBA2_20180701_20181001_MF100_ETOH_r1",
                       "ARMS_Plymouth_MBA2_20180701_20181001_SF40_ETOH_r1",
                       "ARMS_Roscoff_BasBloS1_20180711_20181025_SF40_ETOH_A_r1",
                       "ARMS_Roscoff_BasBloS1_20180711_20181025_SF40_ETOH_B_r1",
                       "ARMS_Roscoff_MarBloR1_20180709_20181024_MF100_ETOH_r1",
                       "ARMS_Roscoff_MarBloR1_20180709_20181024_MF500_ETOH_r1",
                       "ARMS_Roscoff_MarBloR1_20180709_20181024_SF40_ETOH_r1")

ps18S <- prune_samples(!(sample_names(ps18S) %in% to_remove_pres_18S), ps18S)

# Two samples from Roscoff, France have been processed as replicates. #

# Assess rarefaction curves to see which sample to keep
ps_roscoff_18S<-prune_samples(sample_names(ps18S) %in% c("ARMS_Roscoff_BasBloS1_20180711_20181025_SF40_DMSO_A_r1","ARMS_Roscoff_BasBloS1_20180711_20181025_SF40_DMSO_B_r1"), ps18S)
rarecurve(t(otu_table(ps_roscoff_18S)), step=50, cex=0.5)

# Remove the sample with lower sample size/OTU richness
ps18S <- prune_samples(!(sample_names(ps18S) %in% "ARMS_Roscoff_BasBloS1_20180711_20181025_SF40_DMSO_B_r1"), ps18S)

# Remove ASVs which have a total abundance of zero after removing samples during the previous step

ps18S<-prune_taxa(rowSums(otu_table(ps18S))>0,ps18S)

# Check number of OTUs

ps18S

# Get total number of reads

sum(sample_sums(ps18S))

# Get number of OTUs assigned to species level

subset_taxa(ps18S,!is.na(Species))

# Get number of unique species identified with Linnean name 

length(unique(tax_table(subset_taxa(ps18S,!is.na(Species)))[,ncol(tax_table(ps18S))]))

# Get number of species observations for occurrences with a minimum of 2 reads (i.e. all presence-absence occurrences of all OTUs classified to species level across all samples)

ps18S_species<-subset_taxa(ps18S,!is.na(Species))

ps18S_species_obs<-ps18S_species

otu_table(ps18S_species_obs)[otu_table(ps18S_species_obs)<2]<-0

otu_table(ps18S_species_obs)[otu_table(ps18S_species_obs)>1]<-1

sum(sample_sums(ps18S_species_obs))

# get data set agglomerated at phylum level and with relative abundances #

ps18S.phylum <- tax_glom(ps18S, taxrank = "Phylum.Class",NArm=F) 

phylum_18S<-as.data.frame(cbind(as.data.frame(tax_table(ps18S.phylum))[,4],taxa_sums(ps18S.phylum)))

colnames(phylum_18S)<-c("Phylum","reads")

phylum_18S$reads<-as.numeric(phylum_18S$reads)

phylum_18S[is.na(phylum_18S)]<-"NA"

phylum_18S<-aggregate(.~ Phylum,data = phylum_18S,FUN=sum)

phylum_18S<-phylum_18S[order(phylum_18S$reads, decreasing = TRUE),]  

phylum_18S$reads<-phylum_18S$reads/sum(phylum_18S$reads)

# Change classification in the phylum/class level column to correct phylum name
# This was done based on web-based search using WoRMS or other scientific publications

phylum_18S$Phylum # check all unique entries in the Phylum/Class column of the 18S data set

phylum_18S$Phylum<-ifelse(grepl("Urochordata|Craniata|Cephalochordata", phylum_18S$Phylum),"Chordata",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Gregarnimorphea|Dinophyceae|Coccidiomorphea|Syndiniales|Colpodellidea|Perkinsida|Ellobiophyceae|Gregarinomorphea", phylum_18S$Phylum),"Myzozoa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Phaeophyceae|Dictyochophyceae|Chrysophyceae|Synurophyceae|Chrysomerophyceae|Raphidophyceae|Pelagophyceae|Eustigmatophyceae|MOCH-5|Pinguiophyceae|Xanthophyceae", phylum_18S$Phylum),"Ochrophyta",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Florideophyceae|Bangiophyceae|Compsopogonophyceae|Rhodellophyceae|Porphyridiophyceae", phylum_18S$Phylum),"Rhodophyta",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Ulvophyceae|Trebouxiophyceae|Pyramimonadophyceae|Embryophyceae|Chlorophyceae|Mamiellophyceae|Chloropicophyceae|Chlorodendrophyceae|Pedinophyceae|Nephroselmidophyceae|Scotinosphaera", phylum_18S$Phylum),"Chlorophyta",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Spirotrichea|Phyllopharyngea|Oligohymenophorea|Heterotrichea|Litostomatea|CONTH|CONThreeP|Cariacotrichea|Colpodea|Protocruziidae|Prostomatea|Karyorelictea|Nassophorea|Plagiopylea", phylum_18S$Phylum),"Ciliophora",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Labyrinthulomycetes|Bicoecea|Placidideae", phylum_18S$Phylum),"Bigyra",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Developea", phylum_18S$Phylum),"Gyrista",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Thecofilosea|Endomyxa-Ascetosporea|Imbricatea|Endomyxa|Phytomyxea|Filosa-Sarcomonadea|Filosa-Granofilosea|Chlorarachniophyceae|Filosa|Novel-clade-10-12|Metromonadea|Filosa-Thecofilosea|Filosa-Imbricatea", phylum_18S$Phylum),"Cercozoa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Ichthyosporea|Choanoflagellatea", phylum_18S$Phylum),"Choanozoa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Monothalamids|Globothalamea|Allogromida|Tubothalamea", phylum_18S$Phylum),"Foraminifera",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Telonemia-Group-2|Telonemia-Group-1|Katablepharidaceae|Cryptophyceae", phylum_18S$Phylum),"Cryptophyta",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Group-1", phylum_18S$Phylum),"Apusozoa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Myxozoa", phylum_18S$Phylum),"Cnidaria",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Loxomitra", phylum_18S$Phylum),"Entoprocta",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Crasiella", phylum_18S$Phylum),"Gastrotricha",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Prymnesiophyceae|HAP5|Pavlovophyceae", phylum_18S$Phylum),"Haptophyta",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Planomonadida|Subulatomonas-lineage|Mantamonadida|YS16Ec34-lineage", phylum_18S$Phylum),"Sulcozoa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Tubulinea|Discosea-Flabellinia|Stygamoebida|LKM74-lineage|Variosea", phylum_18S$Phylum),"Amoebozoa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Diplonemea|Kinetoplastea|Euglenida|Symbiontida", phylum_18S$Phylum),"Euglenozoa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Prasinodermophyceae", phylum_18S$Phylum),"Prasinodermophyta",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Pterocystida", phylum_18S$Phylum),"Heliozoa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Zygnematophyceae|Zygnemophyceae", phylum_18S$Phylum),"Charophyta",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Pirsonia", phylum_18S$Phylum),"Hyphochytridiomycota",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("RAD-B|Polycystinea", phylum_18S$Phylum),"Radiozoa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Baseodiscus", phylum_18S$Phylum),"Nemertea",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Heterolobosea", phylum_18S$Phylum),"Percolozoa",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Preaxostyla", phylum_18S$Phylum),"Metamonada",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Gnosonesima", phylum_18S$Phylum),"Platyhelminthes",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("twista", phylum_18S$Phylum),"Alveidia",phylum_18S$Phylum)
phylum_18S$Phylum<-ifelse(grepl("Hyphochytridiomycota", phylum_18S$Phylum),"Hyphochytriomyceta",phylum_18S$Phylum)

phylum_18S<-aggregate(.~ Phylum,data = phylum_18S,FUN=sum)

phylum_18S<-phylum_18S[order(phylum_18S$reads, decreasing = TRUE),]

# get number and names of identified phyla
# do not count the MAST clades as unique phyla, they either belong to Bigyra or Gyrista which are already present in the table
phyla_identified_18S<-phylum_18S
length(phyla_identified_18S[!grepl("MAST-",phyla_identified_18S$Phylum),]$Phylum)-1 # Number of recovered phyla; -1 to remove "NA" as phylum
phyla_identified_18S[phyla_identified_18S$Phylum=="NA",] # Percentage of reads unclassified at phylum level
write.table(phyla_identified_18S,"phyla_identified_18S.txt",sep="\t",row.names = F,quote = F)

phylum_18S$Phylum[12:nrow(phylum_18S)]<-"Other" # keep only top 10 most abundant phyla and NA and set rest to "Other"

phylum_18S<-aggregate(.~ Phylum,data = phylum_18S,FUN=sum)

phylum_18S<-phylum_18S[order(phylum_18S$reads, decreasing = TRUE),]

phylum_18S<-phylum_18S[c(1:5,7,9:12,6,8),] # Re-order for plot later on

# For OTUs with species level classification, get data set agglomerated at phylum/class level 

taxa_18S<-as.data.frame(tax_table(ps18S_species))

specs_18S<-distinct(taxa_18S,Species,.keep_all=T)

specs_18S$abundance<-rep(1,nrow(specs_18S))

specs_18S<-specs_18S %>% select(c(Phylum.Class,abundance))

specs_18S<-aggregate(.~ Phylum.Class,data = specs_18S,FUN=sum)

# Change classification in the phylum/class level column to correct phylum name
# This was done based on web-based search using WoRMS or other scientific publications
# We just use the code from above here, because the species-assigned data set is just a subset of the data set above

specs_18S$Phylum.Class<-ifelse(grepl("Urochordata|Craniata|Cephalochordata", specs_18S$Phylum.Class),"Chordata",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Gregarnimorphea|Dinophyceae|Coccidiomorphea|Syndiniales|Colpodellidea|Perkinsida|Ellobiophyceae|Gregarinomorphea", specs_18S$Phylum.Class),"Myzozoa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Phaeophyceae|Dictyochophyceae|Chrysophyceae|Synurophyceae|Chrysomerophyceae|Raphidophyceae|Pelagophyceae|Eustigmatophyceae|MOCH-5|Pinguiophyceae|Xanthophyceae", specs_18S$Phylum.Class),"Ochrophyta",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Florideophyceae|Bangiophyceae|Compsopogonophyceae|Rhodellophyceae|Porphyridiophyceae", specs_18S$Phylum.Class),"Rhodophyta",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Ulvophyceae|Trebouxiophyceae|Pyramimonadophyceae|Embryophyceae|Chlorophyceae|Mamiellophyceae|Chloropicophyceae|Chlorodendrophyceae|Pedinophyceae|Nephroselmidophyceae|Scotinosphaera", specs_18S$Phylum.Class),"Chlorophyta",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Spirotrichea|Phyllopharyngea|Oligohymenophorea|Heterotrichea|Litostomatea|CONTH|CONThreeP|Cariacotrichea|Colpodea|Protocruziidae|Prostomatea|Karyorelictea|Nassophorea|Plagiopylea", specs_18S$Phylum.Class),"Ciliophora",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Labyrinthulomycetes|Bicoecea|Placidideae", specs_18S$Phylum.Class),"Bigyra",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Developea", specs_18S$Phylum.Class),"Gyrista",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Thecofilosea|Endomyxa-Ascetosporea|Imbricatea|Endomyxa|Phytomyxea|Filosa-Sarcomonadea|Filosa-Granofilosea|Chlorarachniophyceae|Filosa|Novel-clade-10-12|Metromonadea", specs_18S$Phylum.Class),"Cercozoa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Ichthyosporea|Choanoflagellatea", specs_18S$Phylum.Class),"Choanozoa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Monothalamids|Globothalamea|Allogromida|Tubothalamea", specs_18S$Phylum.Class),"Foraminifera",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Telonemia-Group-2|Telonemia-Group-1|Katablepharidaceae|Cryptophyceae", specs_18S$Phylum.Class),"Cryptophyta",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Group-1", specs_18S$Phylum.Class),"Apusozoa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Myxozoa", specs_18S$Phylum.Class),"Cnidaria",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Loxomitra", specs_18S$Phylum.Class),"Entoprocta",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Crasiella", specs_18S$Phylum.Class),"Gastrotricha",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Prymnesiophyceae|HAP5|Pavlovophyceae", specs_18S$Phylum.Class),"Haptophyta",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Planomonadida|Subulatomonas-lineage|Mantamonadida|YS16Ec34-lineage", specs_18S$Phylum.Class),"Sulcozoa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Tubulinea|Discosea-Flabellinia|Stygamoebida|LKM74-lineage|Variosea", specs_18S$Phylum.Class),"Amoebozoa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Diplonemea|Kinetoplastea|Euglenida|Symbiontida", specs_18S$Phylum.Class),"Euglenozoa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Prasinodermophyceae", specs_18S$Phylum.Class),"Prasinodermophyta",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Pterocystida", specs_18S$Phylum.Class),"Heliozoa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Zygnematophyceae|Zygnemophyceae", specs_18S$Phylum.Class),"Charophyta",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Pirsonia", specs_18S$Phylum.Class),"Hyphochytridiomycota",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("RAD-B|Polycystinea", specs_18S$Phylum.Class),"Radiozoa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Baseodiscus", specs_18S$Phylum.Class),"Nemertea",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Heterolobosea", specs_18S$Phylum.Class),"Percolozoa",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Preaxostyla", specs_18S$Phylum.Class),"Metamonada",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Gnosonesima", specs_18S$Phylum.Class),"Platyhelminthes",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("twista", specs_18S$Phylum.Class),"Ancoracystidae",specs_18S$Phylum.Class)
specs_18S$Phylum.Class<-ifelse(grepl("Hyphochytridiomycota", specs_18S$Phylum.Class),"Hyphochytriomyceta",specs_18S$Phylum.Class)

specs_18S<-aggregate(.~ Phylum.Class,data = specs_18S,FUN=sum)

specs_18S<-specs_18S[order(specs_18S$abundance, decreasing = TRUE),]

# get number and names of identified phyla
phyla_specs_identified_18S<-specs_18S
length(phyla_specs_identified_18S$Phylum) # Number of recovered phyla
write.table(phyla_specs_identified_18S,"phyla_specs_identified_18S.txt",sep="\t",row.names = F,quote = F)

specs_18S$Phylum.Class<-ifelse(specs_18S$abundance>2,specs_18S$Phylum.Class,"Other") # Keep only phyla with at least 3 species and set rest to "Other"

specs_18S<-aggregate(.~ Phylum.Class,data = specs_18S,FUN=sum)

specs_18S<-specs_18S[order(specs_18S$abundance, decreasing = TRUE),]

specs_18S<-specs_18S[c(1,3:14,2),] # Re-order for plot later on

## ITS ##

# Read ASV counts

ASVcountsITS<-read.table("ITS/ITS_ASV_table.txt",header=T,check.names=F, sep="\t",row.names = 1)

# Read ASV taxonomy (needs to be read as matrix, may cause problems otherwise when phyloseq object will be created)

ASVtaxaITS<-as.matrix(read.table("ITS/ITS_tax_table.txt",header=T,check.names=F, sep="\t",row.names = 1))

# Sort count table based on order in tax table (precautionary measure)

ASVcountsITS<-ASVcountsITS[order(match(rownames(ASVcountsITS),rownames(ASVtaxaITS))),]

# Create phyloseq object

psITS <- phyloseq(otu_table(ASVcountsITS,taxa_are_rows = TRUE), sample_data(samples), tax_table(ASVtaxaITS))

# Get number of samples that produced ASVs through PEMA processing

psITS

# Remove the sediment and plankton samples (some sediment and plankton samples were sequenced as a trial during the initial phase of the ARMS program)

psITS <- subset_samples(psITS,Fraction!="SED")
psITS <- subset_samples(psITS,Fraction!="PS")

# One ASV got assigned to a fungal species, but its taxonomy seemed to be incorrectly deposited in the Unite database as a insect genus of the same name
# Replace insect taxonomy with the correct fungal taxonomy of this species

tax_table(psITS)[which(tax_table(psITS)[,9]=="Petrophila_incerta"),3:8]<-c("Fungi","Ascomycota","Dothideomycetes","Mycosphaerellales","Extremaceae","Petrophila")

# Some samples have been preserved in DMSO as well as EtOH initially as a trial. 
# Where for the remaining samples both replicates exist, keep only DMSO samples.

sort(sample_names(psITS)) # Check sample names to see where both replicates exist (check for same sample name with and without "_EtOH"

to_remove_pres_ITS <-c("ARMS_Crete_1HERP_20180928_20190128_MF500_ETOH_r1",
                       "ARMS_Crete_1HERP_20180928_20190128_SF_ETOH_r1",
                       "ARMS_Koster_VH2_20180418_20180906_SF40_ETOH_r1",
                       "ARMS_Plymouth_MBA2_20180701_20181001_MF100_ETOH_r1")

psITS <- prune_samples(!(sample_names(psITS) %in% to_remove_pres_ITS), psITS)

# Remove ASVs which have a total abundance of zero after removing samples during the previous step

psITS<-prune_taxa(rowSums(otu_table(psITS))>0,psITS)

# Check number of ASVs

psITS

# Get total number of reads

sum(sample_sums(psITS))

# Get number of ASVs assigned to species level

subset_taxa(psITS,!is.na(Species))

# Get number of unique species identified with Linnean name 

length(unique(tax_table(subset_taxa(psITS,!is.na(Species)))[,ncol(tax_table(psITS))]))

# Get number of species observations for occurrences with a minimum of 2 reads (i.e. all presence-absence occurrences of all ASVs classified to species level across all samples)

psITS_species<-subset_taxa(psITS,!is.na(Species))

psITS_species_obs<-psITS_species

otu_table(psITS_species_obs)[otu_table(psITS_species_obs)<2]<-0

otu_table(psITS_species_obs)[otu_table(psITS_species_obs)>1]<-1

sum(sample_sums(psITS_species_obs))

# get data set agglomerated at phylum level and with relative abundances

psITS.phylum <- tax_glom(psITS, taxrank = "Phylum",NArm=F) 

phylum_ITS<-as.data.frame(cbind(as.data.frame(tax_table(psITS.phylum))[,4],taxa_sums(psITS.phylum)))

colnames(phylum_ITS)<-c("Phylum","reads")

phylum_ITS$reads<-as.numeric(phylum_ITS$reads)

phylum_ITS[is.na(phylum_ITS)]<-"NA"

phylum_ITS<-aggregate(.~ Phylum,data = phylum_ITS,FUN=sum)

phylum_ITS<-phylum_ITS[order(phylum_ITS$reads, decreasing = TRUE),]  

phylum_ITS$reads<-phylum_ITS$reads/sum(phylum_ITS$reads)

# All entries in the Phylum level column represent correct phylum names. No correction needed.
phylum_ITS$Phylum 

# get number and names of identified phyla
phyla_identified_ITS<-phylum_ITS
length(which(phyla_identified_ITS$Phylum!="NA")) # Get number of phyla 
phyla_identified_ITS[phyla_identified_ITS$Phylum=="NA",] # Percentage of reads unclassified at phylum level
write.table(phyla_identified_ITS,"phyla_identified_ITS.txt",sep="\t",row.names = F,quote = F)

phylum_ITS$Phylum[7:nrow(phylum_ITS)]<-"Other" # Keep only top 5 most abundant phyla and set rest to "Other"

phylum_ITS<-aggregate(.~ Phylum,data = phylum_ITS,FUN=sum)

phylum_ITS<-phylum_ITS[order(phylum_ITS$reads, decreasing = TRUE),]

phylum_ITS<-phylum_ITS[c(1,3:6,7,2),] # Re-order for plot later on

# For ASVs with species level classification, get data set agglomerated at class  level.
# we chose class level here for a better resolution in the plot.

taxa_ITS<-as.data.frame(tax_table(psITS_species))

specs_ITS<-distinct(taxa_ITS,Species,.keep_all=T)

specs_ITS$abundance<-rep(1,nrow(specs_ITS))

specs_ITS<-specs_ITS %>% select(c(Class,abundance))

specs_ITS<-aggregate(.~ Class,data = specs_ITS,FUN=sum)

specs_ITS<-specs_ITS[order(specs_ITS$abundance, decreasing = TRUE),]

# All entries in the class level column represent correct class names. No correction needed.

# get number and names of identified classes and manually add phylum classification for each class
phyla_specs_identified_ITS<-specs_ITS
length(phyla_specs_identified_ITS$Class) # Number of recovered classes
phyla_specs_identified_ITS$Phylum<-"unknown" # Prepare to add respective phylum classification for the classes
phyla_specs_identified_ITS$Phylum<-ifelse(grepl("Eurotiomycetes|Dothideomycetes|Sordariomycetes|Saccharomycetes|Leotiomycetes|Lecanoromycetes|Pezizomycetes", phyla_specs_identified_ITS$Class),"Ascomycota",phyla_specs_identified_ITS$Phylum)
phyla_specs_identified_ITS$Phylum<-ifelse(grepl("Agaricomycetes|Tremellomycetes|Agaricostilbomycetes|Microbotryomycetes|Wallemiomycetes|Malasseziomycetes", phyla_specs_identified_ITS$Class),"Basidiomycota",phyla_specs_identified_ITS$Phylum)
phyla_specs_identified_ITS$Phylum<-ifelse(grepl("Mortierellomycetes", phyla_specs_identified_ITS$Class),"Mortierellomycota",phyla_specs_identified_ITS$Phylum)
phyla_specs_identified_ITS$Phylum<-ifelse(grepl("Mucoromycetes", phyla_specs_identified_ITS$Class),"Mucoromycota",phyla_specs_identified_ITS$Phylum)
write.table(phyla_specs_identified_ITS,"phyla_class_specs_identified_ITS.txt",sep="\t",row.names = F,quote = F)

specs_ITS$Class<-ifelse(specs_ITS$abundance>2,specs_ITS$Class,"Other") # Keep only phyla with at least 3 species and set rest to "Other"

specs_ITS<-aggregate(.~ Class,data = specs_ITS,FUN=sum)

specs_ITS<-specs_ITS[order(specs_ITS$abundance, decreasing = TRUE),]

specs_ITS<-specs_ITS[c(1:3,5:8,4),] # Re-order for plot later on

### Get overall info on phyla and species recovered ###

phyla<-unique(c(phyla_identified_COI$Phylum,phyla_identified_18S$Phylum,phyla_identified_ITS$Phylum)) # vector with all unique phyla of all genes
length(which(!grepl("MAST-",phyla)))-1 # Number of recovered phyla excluding the MASTs (they are already included in either Bigyra or Gyrista phyla); -1 to remove "NA" as phylum

species<-unique(c(phyla_specs_identified_COI$Phylum,phyla_specs_identified_18S$Phylum,phyla_specs_identified_ITS$Phylum)) # vector with all unique species of all genes
length(species)

### Plots of phyla recovered and phylum/class classification of species recovered ###

# Create table with colors for each Phylum/Class 

phylum_palette<-as.data.frame(c(phylum_COI[,1],phylum_18S[,1],phylum_ITS[,1],specs_ITS[,1],specs_18S[,1],specs_COI[,1]))

phylum_palette<-unique(phylum_palette)

write.table(phylum_palette,"palette_to_be_formatted.txt",sep = "\t",col.names = F,row.names = F,quote=F)

# Check kelly palette of grafify package for colorblind-friendly codes

graf_palettes$kelly

# outside of R, add color codes as a second column to palette_to_be_formatted.txt in Excel. Give the fungal classes colors according to the phylum they belong to
# the kelly palette contains less colors than needed, choose a few more suitable hexacodes and add them manually
# Save file as palette.txt once ready
# try out plotting with the colors and switch hexacodes around as needed for a nice order of colors in the plots

# read prepared palette table back into R

phylum_palette<-read.table("palette.txt",sep="\t",stringsAsFactors = FALSE,na.strings = "",comment.char = "")

colnames(phylum_palette)<-c("Phylum","color")

# Plot with gene names on the left side

cols1<-phylum_palette[phylum_palette$Phylum %in% phylum_COI$Phylum,]
cols1<-cols1[order(match(cols1[,1],phylum_COI$Phylum)),]
cols1<-cols1$color
phylum_COI_plot<-ggplot(phylum_COI, aes(x=Phylum, y=reads, fill=Phylum))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=phylum_COI$Phylum)+
  scale_y_continuous(limits=c(0,0.5),breaks = c(0,.1,.2,.3,.4,.5))+
  scale_fill_manual(breaks=phylum_COI$Phylum,values=cols1)+
  ggtitle("A")+
  annotate("text", y = 0.25, x = -4, label = "COI", size = 7, fontface =2)+
  coord_cartesian(xlim = c(1,12),clip="off")+
  theme_bw()+
  theme(plot.margin=unit(c(.2,.5,0,2.2), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",axis.text.y=element_text(color="black",size=10),axis.text.x=element_text(angle=45,hjust=1,color="black",size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank())

cols2<-phylum_palette[phylum_palette$Phylum %in% phylum_18S$Phylum,]
cols2<-cols2[order(match(cols2[,1],phylum_18S$Phylum)),]
cols2<-cols2$color
phylum_18S_plot<-ggplot(phylum_18S, aes(x=Phylum, y=reads, fill=Phylum))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=phylum_18S$Phylum)+
  scale_y_continuous(limits=c(0,0.3),breaks = c(0,.1,.2,.3))+
  scale_fill_manual(breaks=phylum_18S$Phylum,values=cols2)+
  ggtitle("C")+
  annotate("text", y = 0.15, x = -4, label = "18S", size = 7, fontface =2)+
  coord_cartesian(xlim = c(1,12),clip="off")+
  labs(y="Relative abundance")+
  theme_bw()+
  theme(plot.margin=unit(c(-0.3,.5,-0.3,2.2), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.x=element_blank(),legend.position="none",axis.text.y=element_text(color="black",size=10),axis.text.x=element_text(angle=45,hjust=1,color="black",size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank(),axis.title.y=element_text(size=18))

cols3<-phylum_palette[phylum_palette$Phylum %in% phylum_ITS$Phylum,]
cols3<-cols3[order(match(cols3[,1],phylum_ITS$Phylum)),]
cols3<-cols3$color
phylum_ITS_plot<-ggplot(phylum_ITS, aes(x=Phylum, y=reads, fill=Phylum))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=phylum_ITS$Phylum)+
  scale_y_continuous(limits=c(0,0.45),breaks = c(0,0.2,0.4))+
  scale_fill_manual(breaks=phylum_ITS$Phylum,values=cols3)+
  ggtitle("E")+
  annotate("text", y = 0.3, x = -2.2, label = "ITS", size = 7, fontface =2)+
  coord_cartesian(xlim = c(1,7),clip="off")+
  theme_bw()+
  theme(plot.margin=unit(c(0,.5,.2,2.2), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position="none",axis.text.y=element_text(color="black",size=10),axis.text.x=element_text(angle=45,hjust=1,color="black",size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank())

cols4<-phylum_palette[phylum_palette$Phylum %in% specs_COI$Phylum,]
cols4<-cols4[order(match(cols4[,1],specs_COI$Phylum)),]
cols4<-cols4$color
specs_COI_plot<-ggplot(specs_COI, aes(x=Phylum, y=abundance, fill=Phylum))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=specs_COI$Phylum)+
  scale_y_continuous(limits=c(0,130),breaks = c(0,30,60,90,120))+
  scale_fill_manual(breaks=specs_COI$Phylum,values=cols4)+
  ggtitle("B")+
  theme_bw()+
  theme(plot.margin=unit(c(.2,.2,0,0), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",axis.text.y=element_text(color="black",size=10),axis.text.x=element_text(angle=45,hjust=1,color="black",size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank())

cols5<-phylum_palette[phylum_palette$Phylum %in% specs_18S$Phylum,]
cols5<-cols5[order(match(cols5[,1],specs_18S$Phylum)),]
cols5<-cols5$color
specs_18S_plot<-ggplot(specs_18S, aes(x=Phylum.Class, y=abundance, fill=Phylum.Class))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=specs_18S$Phylum.Class)+
  scale_y_continuous(limits=c(0,40),breaks = c(0,10,20,30,40))+
  labs(y="Number of species")+
  ggtitle("D")+
  scale_fill_manual(breaks=specs_18S$Phylum,values=cols5)+
  theme_bw()+
  theme(plot.margin=unit(c(-0.3,.2,-0.3,0), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.x=element_blank(),legend.position="none",axis.text.y=element_text(color="black",size=10),axis.text.x=element_text(angle=45,hjust=1,color="black",size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank(),axis.title.y=element_text(size=18))

cols6<-phylum_palette[phylum_palette$Phylum %in% specs_ITS$Class,]
cols6<-cols6[order(match(cols6[,1],specs_ITS$Class)),]
cols6<-cols6$color
specs_ITS_plot<-ggplot(specs_ITS, aes(x=Class, y=abundance, fill=Class))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=specs_ITS$Class)+
  scale_fill_manual(breaks=specs_ITS$Class,values=cols6)+
  ggtitle("F")+
  theme_bw()+
  theme(plot.margin=unit(c(0,.2,.2,0), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",axis.text.y=element_text(color="black",size=10),axis.text.x=element_text(angle=45,hjust=1,color="black",size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank())

phylum_plot_genes<-egg::ggarrange(phylum_COI_plot,specs_COI_plot,phylum_18S_plot,specs_18S_plot,phylum_ITS_plot,specs_ITS_plot,nrow=3)

ggsave(phylum_plot_genes,file="Figures/Figure_2.pdf",height=8,width=8,device = cairo_pdf,dpi=600)

# Plot without gene names

cols1<-phylum_palette[phylum_palette$Phylum %in% phylum_COI$Phylum,]
cols1<-cols1[order(match(cols1[,1],phylum_COI$Phylum)),]
cols1<-cols1$color
phylum_COI_plot<-ggplot(phylum_COI, aes(x=Phylum, y=reads, fill=Phylum))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=phylum_COI$Phylum)+
  scale_y_continuous(limits=c(0,0.5),breaks = c(0,.1,.2,.3,.4,.5))+
  scale_fill_manual(breaks=phylum_COI$Phylum,values=cols1)+
  ggtitle("A")+
  theme_bw()+
  theme(plot.margin=unit(c(0,0,-0.3,0), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",axis.text.y=element_text(color="black"),axis.text.x=element_text(angle=45,hjust=1,color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank())

cols2<-phylum_palette[phylum_palette$Phylum %in% phylum_18S$Phylum,]
cols2<-cols2[order(match(cols2[,1],phylum_18S$Phylum)),]
cols2<-cols2$color
phylum_18S_plot<-ggplot(phylum_18S, aes(x=Phylum, y=reads, fill=Phylum))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=phylum_18S$Phylum)+
  scale_y_continuous(limits=c(0,0.3),breaks = c(0,.1,.2,.3))+
  scale_fill_manual(breaks=phylum_18S$Phylum,values=cols2)+
  ggtitle("C")+
  labs(y="Relative abundance")+
  theme_bw()+
  theme(plot.margin=unit(c(0,0,-0.3,0), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.x=element_blank(),legend.position="none",axis.text.y=element_text(color="black"),axis.text.x=element_text(angle=45,hjust=1,color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank(),axis.title.y=element_text(size=15))

cols3<-phylum_palette[phylum_palette$Phylum %in% phylum_ITS$Phylum,]
cols3<-cols3[order(match(cols3[,1],phylum_ITS$Phylum)),]
cols3<-cols3$color
phylum_ITS_plot<-ggplot(phylum_ITS, aes(x=Phylum, y=reads, fill=Phylum))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=phylum_ITS$Phylum)+
  scale_y_continuous(limits=c(0,0.45),breaks = c(0,0.2,0.4))+
  scale_fill_manual(breaks=phylum_ITS$Phylum,values=cols3)+
  ggtitle("E")+
  theme_bw()+
  theme(plot.margin=unit(c(0,0,0,0), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position="none",axis.text.y=element_text(color="black"),axis.text.x=element_text(angle=45,hjust=1,color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank())

cols4<-phylum_palette[phylum_palette$Phylum %in% specs_COI$Phylum,]
cols4<-cols4[order(match(cols4[,1],specs_COI$Phylum)),]
cols4<-cols4$color
specs_COI_plot<-ggplot(specs_COI, aes(x=Phylum, y=abundance, fill=Phylum))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=specs_COI$Phylum)+
  scale_y_continuous(limits=c(0,130),breaks = c(0,30,60,90,120))+
  scale_fill_manual(breaks=specs_COI$Phylum,values=cols4)+
  ggtitle("B")+
  theme_bw()+
  theme(plot.margin=unit(c(0,0,-0.3,.5), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",axis.text.y=element_text(color="black"),axis.text.x=element_text(angle=45,hjust=1,color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank())

cols5<-phylum_palette[phylum_palette$Phylum %in% specs_18S$Phylum,]
cols5<-cols5[order(match(cols5[,1],specs_18S$Phylum)),]
cols5<-cols5$color
specs_18S_plot<-ggplot(specs_18S, aes(x=Phylum.Class, y=abundance, fill=Phylum.Class))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=specs_18S$Phylum.Class)+
  scale_y_continuous(limits=c(0,40),breaks = c(0,10,20,30,40))+
  labs(y="Number of species")+
  ggtitle("D")+
  scale_fill_manual(breaks=specs_18S$Phylum,values=cols5)+
  theme_bw()+
  theme(plot.margin=unit(c(0,0,-0.3,.5), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.x=element_blank(),legend.position="none",axis.text.y=element_text(color="black"),axis.text.x=element_text(angle=45,hjust=1,color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank(),axis.title.y=element_text(size=15))

cols6<-phylum_palette[phylum_palette$Phylum %in% specs_ITS$Class,]
cols6<-cols6[order(match(cols6[,1],specs_ITS$Class)),]
cols6<-cols6$color
specs_ITS_plot<-ggplot(specs_ITS, aes(x=Class, y=abundance, fill=Class))+
  geom_bar(stat="identity", color="black")+
  scale_x_discrete(limits=specs_ITS$Class)+
  scale_fill_manual(breaks=specs_ITS$Class,values=cols6)+
  ggtitle("F")+
  theme_bw()+
  theme(plot.margin=unit(c(0,0,0,.5), "cm"),plot.title = element_text(size = 20, face = "bold",vjust=-.2),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none",axis.text.y=element_text(color="black"),axis.text.x=element_text(angle=45,hjust=1,color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),panel.border = element_blank())

phylum_plot<-egg::ggarrange(phylum_COI_plot,specs_COI_plot,phylum_18S_plot,specs_18S_plot,phylum_ITS_plot,specs_ITS_plot,nrow=3)

ggsave(phylum_plot,file="phylum_plot.png",height=8,width=8,type="cairo")

### Get species occurrences to scan species lists against IUCN/HELCOm red lists, AMBI and WRiMS ###

# Make list of species to get accepted name and AphiaID using LifeWatch Belgium's E-services 

species_ps_COI<-tax_glom(psCOI_species_obs,"Species")
species_ps_COI<-merge_samples(species_ps_COI,"ARMS")
otu_table(species_ps_COI)<-t(otu_table(species_ps_COI))
species_list_COI_all<-as.data.frame(cbind(otu_table(species_ps_COI),tax_table(species_ps_COI)[,7]))
species_list_COI_all<-species_list_COI_all %>% pivot_longer(cols=-Species,names_to = "ARMS", values_to = "count")

species_ps_18S<-tax_glom(ps18S_species_obs,"Species")
species_ps_18S<-merge_samples(species_ps_18S,"ARMS")
otu_table(species_ps_18S)<-t(otu_table(species_ps_18S))
species_list_18S_all<-as.data.frame(cbind(otu_table(species_ps_18S),tax_table(species_ps_18S)[,10]))
species_list_18S_all<-species_list_18S_all %>% pivot_longer(cols=-Species,names_to = "ARMS", values_to = "count")

species_list<-rbind(species_list_COI_all,species_list_18S_all)
species_list<-as.data.frame(unique(species_list$Species))
colnames(species_list)<-"ScientificName" # The LifeWatch E-servcies used later on needs the Species column to be named "ScientificName"
species_list$ScientificName<-gsub("_"," ",species_list$ScientificName)
write.table(species_list,"species_list.txt",row.names = F,sep="\t",quote=F)

# Go to LifeWatch Belgium's e-Lab services: https://www.lifewatch.be/data-services/
# Under "Run Services", use the tables just created as input  and run "Taxon match services" -> "Taxon match World Register of Marine Species (WoRMS)"
# This will give us AphiaIDs for the species which are needed to scan against the databases
# The following steps were done outside of R in Excel:
# For cases where there where so-called "fuzzy matches" and therefore no AphiaID was returned, we manually looked up the species in WoRMS, chose the most matching synonymised entry based on distribution, and entered the respective values in the columns accepted_name_aphia_worms and valid_aphiaid_worms.
# For cases where WoRMS status was "uncertain", we manually added the taxon name AphiaID present in WoRMS to the column accepted_name_aphia_worms and valid_aphiaid_worms.
# For cases were there was no match in WoRMS, the taxon was removed.
# The table was then subset to the following columns: scientificname,accepted_name_aphia_worms and valid_aphiaid_worms

# Read the prepared results file from WoRMS Taxon Match procedure run on LifeWatch E-services 
# Add accepted name and AphiaID to the COI and 18S species lists created above and safe files used as input for AMBI and IUCN/HELCOM Red List screening

results_species<-read.table("result_species_list.txt",header=T,sep="\t")

species_list_COI_all$count<-as.numeric(species_list_COI_all$count)
species_list_18S_all$count<-as.numeric(species_list_18S_all$count)

species_list_COI<-species_list_COI_all[species_list_COI_all$count>0,]
species_list_18S<-species_list_18S_all[species_list_18S_all$count>0,]

species_list_COI$Species<-gsub("_"," ",species_list_COI$Species)
species_list_18S$Species<-gsub("_"," ",species_list_18S$Species)

species_list_COI<-species_list_COI %>% filter(Species %in% results_species$scientificname)
species_list_18S<-species_list_18S %>% filter(Species %in% results_species$scientificname)

for(i in 1:nrow(results_species)) { # add accepted species name and AphiID
  species_list_COI$accepted_name_aphia_worms[species_list_COI$Species == results_species$scientificname[i]] <- results_species$accepted_name_aphia_worms[i]
  species_list_COI$valid_aphiaid_worms[species_list_COI$accepted_name_aphia_worms == results_species$accepted_name_aphia_worms[i]] <- results_species$valid_aphiaid_worms[i]
  species_list_18S$accepted_name_aphia_worms[species_list_18S$Species == results_species$scientificname[i]] <- results_species$accepted_name_aphia_worms[i]
  species_list_18S$valid_aphiaid_worms[species_list_18S$accepted_name_aphia_worms == results_species$accepted_name_aphia_worms[i]] <- results_species$valid_aphiaid_worms[i]
} # will give a warning, just ignore

for(i in 1:nrow(samples)) { # add Observatory info
  species_list_COI$Observatory[species_list_COI$ARMS == samples$ARMS[i]] <- samples$Observatory[i]
  species_list_18S$Observatory[species_list_18S$ARMS == samples$ARMS[i]] <- samples$Observatory[i]
} # will give a warning, just ignore

species_list_COI<-species_list_COI %>% select(-c(Species,count)) %>% relocate(accepted_name_aphia_worms,valid_aphiaid_worms)
species_list_18S<-species_list_18S %>% select(-c(Species,count)) %>% relocate(accepted_name_aphia_worms,valid_aphiaid_worms)

colnames(species_list_COI)[1:2]<-c("Species","AphiaID")
colnames(species_list_18S)[1:2]<-c("Species","AphiaID")

species_list_COI<-species_list_COI[order(species_list_COI$ARMS),]
species_list_18S<-species_list_18S[order(species_list_18S$ARMS),]

write.csv(species_list_COI,"SpeciesListCOI.csv",row.names = F)
write.csv(species_list_18S,"SpeciesList18S.csv",row.names = F)
# IUCN/HELCOM Red List and AMBI scan was then performed with these .csv files #

# Make table for WRiMS scan #

species_list_wrims<-rbind(species_list_COI,species_list_18S)
species_list_wrims$count<-1 # add column with occurrence cout of 1 for presence
species_list_wrims<-species_list_wrims %>% distinct()

# Create list with separate tables for each ARMS

species_list_wrims <- species_list_wrims %>%  group_by(ARMS) # make grouped data frame
group_name_df <- group_keys(species_list_wrims) %>%  mutate(group_name = ARMS) # get group keys
group_name <- group_name_df$group_name # get name for each group
species_list_wrims <- group_split(species_list_wrims) %>% setNames(group_name) # assign name to each split table
species_list_wrims<-lapply(species_list_wrims,function(x) {colnames(x)[5]<-x[1,3];x}) # Rename count with name of respective ARMS unit
species_list_wrims<-lapply(species_list_wrims,function(x) x %>% select(-c(ARMS,Observatory))) # remove columns

species_list_wrims<-lapply(species_list_wrims, function(x) {x$NCBI_Tax.x_pr2<-"dummy:12345";x}) # Add "NCBI_Tax.x_pr2" column with dummy values, is needed for the Jupyter notebook code
species_list_wrims<-lapply(species_list_wrims, function(x) {x$ScientificName_pr2<-"dummy dummy";x})  # Add "ScientificName_pr2" column with dummy values, is needed for the Jupyter notebook code
species_list_wrims<-lapply(species_list_wrims, function(x) x %>% relocate(c(NCBI_Tax.x_pr2,ScientificName_pr2,AphiaID)))
species_list_wrims<-lapply(species_list_wrims, function(x) {colnames(x)[3:4]<-c("AphiaID_accepted","ScientificName_accepted");x}) # rename columns as required for Junyper notebook code
species_list_wrims<-lapply(species_list_wrims, as.data.frame) # as.data.frame necessary to avoid error when using xlsx::write.xlsx with row.names=F

# Write xlsx file with all list elements of species_list_wrims as separate sheets in this file
# Note that the "_18S" ending in the file names was added because the code used for WRiMS scan was initially written requiring this exact file name. The file created here contains both COI AND 18S species occurrences.  
openxlsx::write.xlsx(species_list_wrims, file = "ARMS_SpeciesPerObservatory_18S.xlsx",sheetName = names(species_list_wrims), rowNames = FALSE)

# read Results from AMBI and IUCN scan

iucn_ambi_coi<-read.csv("SpeciesListAttributesCOI.csv",header=T,strip.white=TRUE) # strip.white necessary to remove white spaces
iucn_ambi_18s<-read.csv("SpeciesListAttributes18S.csv",header=T,strip.white=TRUE)

# Subset species listed in: IUCN Red List and HELCOM Red List with known status other than Least concern or Data Deficient, and in AMBI as "very sensitive to disturbance" 
# merge COI and 18S data
# per observatory and anthropogenic influence category, count number of species which are present in these databases

iucn_coi<-iucn_ambi_coi %>% filter(!IUCN.RedList.Category %in% c("None","Data Deficient","Least Concern"))
helcom_coi<-iucn_ambi_coi %>% filter(!HELCOM.RedList.Category %in% c("None","Data Deficient","Least Concern"))
iucn_18s<-iucn_ambi_18s %>% filter(!IUCN.RedList.Category %in% c("None","Data Deficient","Least Concern"))
helcom_18s<-iucn_ambi_18s %>% filter(!HELCOM.RedList.Category %in% c("None","Data Deficient","Least Concern"))
red_list<-rbind(iucn_coi,helcom_coi,iucn_18s,helcom_18s)
red_list<-red_list[red_list$Species.name!="Pinna nobilis",] # remove Pinna nobilis, cause it occurs with 5 reads in Plymouth, but is not known from outside Mediterranean
red_list_n<-red_list %>% group_by(Observatory) %>% summarise(n_red_list=n_distinct(Species.name))

ambi_coi<-iucn_ambi_coi %>% filter(AMBI.ecological.group %in% "very sensitive to disturbance")
ambi_18s<-iucn_ambi_18s %>% filter(AMBI.ecological.group %in% "very sensitive to disturbance")
ambi_list<-rbind(ambi_coi,ambi_18s)
ambi_list_n<-ambi_list %>% group_by(Observatory) %>% summarise(n_ambi_list=n_distinct(Species.name))

# Read results from WRiMS scan as list of all ARMS sheets 

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
wrims_sheets <- read_excel_allsheets("ARMS_SpeciesPerObservatory_wrims.xlsx")

wrims_sheets <- lapply(wrims_sheets, function(x) {colnames(x)[6]<-"occurrence";x}) # column with presence-absence is named after each individual ARMS unit, needs to be changed before combining all list elements
wrims_sheets$MarBloR1<-wrims_sheets$MarBloR1 %>% select(-c(MarBloR2,MarBloR3)) # For some reason, the MarBlo dataframe has two extra columns for the other two MarBlo ARMS: remove them
wrims<-bind_rows(wrims_sheets, .id = "ARMS") # Combine all dataframes into one single dataframe

for(j in 1:nrow(samples)) { # add observatory info for each ARMS
    wrims$Observatory[wrims$ARMS == samples$ARMS[j]] <- samples[,"Observatory"][j]
}

wrims<-wrims[wrims$occurrence!=0,]
wrims<-wrims[which(wrims$establishmentMeans=="Alien"),]
wrims<-wrims %>% select(c(ScientificName_accepted,Observatory)) %>% distinct()
wrims_n<-wrims %>% group_by(Observatory) %>% summarise(n_wrims=n_distinct(ScientificName_accepted))

# Save info on occurring species found in each database at each observatory to file 

ambi_list$list<-"AMBI"
colnames(ambi_list)[1]<-"Species"
length(unique(ambi_list$Species)) # Number of unique species

red_list$list<-"Red_List"
colnames(red_list)[1]<-"Species"
length(unique(red_list$Species)) # Number of unique species

wrims$list<-"WRiMS"
colnames(wrims)[1]<-"Species"
length(unique(wrims$Species)) # Number of unique species

all_list<-rbind(select(ambi_list,c(Species,Observatory,list)),select(wrims,c(Species,Observatory,list)),select(red_list,c(Species,Observatory,list)))

all_list<-all_list %>% distinct()
              
all_list<-all_list %>% group_by(list) %>% # First have to create a unique identifier row for each list type before using pivot_wider
                       mutate(row = row_number()) %>%
                       pivot_wider(names_from = list, values_from = Species) %>%
                       select(-row)
all_list<-all_list %>%
  group_by(Observatory) %>%
  mutate(across(.cols = WRiMS, ~WRiMS[order(is.na(.))])) %>% # This shifts vales in WRiMS column up for each Observatory so NAs are at bottom
  mutate(across(.cols = Red_List, ~Red_List[order(is.na(.))])) # This shifts vales in Red List column up for each Observatory so NAs are at bottom

writexl::write_xlsx(split(subset(all_list, select = -Observatory), all_list$Observatory), "observatories_AMBI_RedList_WRiMS_species.xlsx")

# merge AMBI, Red List and WRiMS tables

ambi_red_list_wrims<-merge(red_list_n,ambi_list_n,all=T) # can only merge two elements at once
ambi_red_list_wrims<-merge(ambi_red_list_wrims,wrims_n,all=T)
ambi_red_list_wrims[is.na(ambi_red_list_wrims)]<-0
ambi_red_list_wrims$arms<-c(2,2,5,4,2,3,1,6,3,5,5,8,4,2,4) # Manually add number of ARMS units per Observatory to order x-axis of plot later on
ambi_red_list_wrims<-ambi_red_list_wrims[order(ambi_red_list_wrims$arms),] # Order by number of ARMs units for plot
write.table(select(ambi_red_list_wrims,-arms),"observatories_AMBI_RedList_WRiMS_counts.txt",sep="\t",row.names=F)
ambi_red_list_wrims<-ambi_red_list_wrims %>% pivot_longer(-c(Observatory,arms))

# Plots

graf_palettes$okabe_ito # colors
ambi_red_list_wrims$name<-factor(ambi_red_list_wrims$name,levels=c("n_ambi_list","n_wrims","n_red_list")) # Re-order for plot

# Plot without legend
species_plot<-ggplot(ambi_red_list_wrims, aes(x = Observatory,y=value,fill=name)) +
  geom_bar(position="dodge", stat="identity",color="black")+
  theme_classic()+
  scale_y_continuous(limits = c(0,40),breaks=c(0,10,20,30,40))+
  scale_x_discrete(expand=expansion(mult = c(0, 0),add = c(2,0)),limits=ambi_red_list_wrims$Observatory)+
  scale_fill_manual(values=c("#56B4E9","#009E73","#D55E00"),labels=c("Species registered as very sensitive to disturbance in AMBI","Species registered as alien at the location of occurrence in WRiMS","Species registered as Near Threatened, Vulnerable, Endangered or Critically Endangered in IUCN and/or HELCOM Red List"))+
  theme(legend.position = "none",axis.title.x = element_blank(),axis.text.y=element_text(color="black"),panel.background = element_rect(fill = "transparent",colour = NA_character_),plot.background = element_rect(fill = "transparent",colour = NA_character_),axis.text.x=element_text(angle=45,hjust=1,color="black"))+
  labs(y="Number of species")
ggsave(species_plot,file="Figures/Figure_6.png",height=3,width=5,dpi=600)

# Plot with legend just to get legend and add it to previous plot outside of R
species_plot_legend<-ggplot(ambi_red_list_wrims, aes(x = Observatory,y=value,fill=name)) +
  geom_bar(position="dodge", stat="identity",color="black")+
  theme_classic()+
  scale_y_continuous(limits = c(0,40),breaks=c(0,10,20,30,40))+
  scale_x_discrete(expand=expansion(mult = c(0, 0),add = c(2,0)),limits=ambi_red_list_wrims$Observatory)+
  scale_fill_manual(values=c("#56B4E9","#009E73","#D55E00"),labels=c("Species registered as very sensitive to disturbance in AMBI","Species registered as alien at the location of occurrence in WRiMS","Species registered as Near Threatened, Vulnerable, Endangered or\nCritically Endangered in IUCN and/or HELCOM Red List"))+
  theme(legend.background = element_rect(fill='transparent'),legend.spacing.y = unit(.3, 'cm'),legend.key.size=unit(.3, "cm"),legend.text = element_text(size=8),legend.position = "bottom",legend.title=element_blank(),axis.title.x = element_blank(),axis.text.y=element_text(color="black"),panel.background = element_rect(fill = "transparent",colour = NA_character_),plot.background = element_rect(fill = "transparent",colour = NA_character_),axis.text.x=element_text(angle=45,hjust=1,color="black"))+
  labs(y="Number of species",fill="Trait category")+
  guides(fill=guide_legend(byrow = TRUE,ncol=1))
ggsave(species_plot_legend,file="Figures/Figure_6_legend.png",height=4,width=6,dpi=600)

# Outside of R, using any image processing software, manually change legend key size for Red list key, change space between legend items and add legend to first plot, then save as pdf.

### UpSet plot for number of species identified which are shared between the marker gene data sets ###

# Get unique species occurrences for each gene

list_coi<-as.data.frame(unique(as.data.frame(tax_table(psCOI_species))$Species))
list_coi$COI<-1
colnames(list_coi)[1]<-"Species"

list_18s<-as.data.frame(unique(as.data.frame(tax_table(ps18S_species))$Species))
list_18s$'18S'<-1
colnames(list_18s)[1]<-"Species"

list_its<-as.data.frame(unique(as.data.frame(tax_table(psITS_species))$Species))
list_its$ITS<-1
colnames(list_its)[1]<-"Species"

# merge tables by species name, set NAs to zero and set species names as rownames

list_all<-merge(list_coi,list_18s,by="Species",all=T)
list_all<-merge(list_all,list_its,by="Species",all=T)

list_all[is.na(list_all)] <- 0
rownames(list_all) <- list_all$Species
list_all$Species <- NULL

# Create UpSet plot and save it to file

png(file="upset.png",height = 1500, width = 1200,res=300,type="cairo")
upset(list_all, main.bar.color = "#009E73", sets.bar.color = "lightgray", matrix.color = "black", order.by = "freq", decreasing = T, set_size.show = T,set_size.scale_max=890,
      mainbar.y.label = "No. of species", sets.x.label = "No. of species",text.scale = c(1.3, 1.3, 1, 1, 1.3, 1.5))
dev.off()

### Sequencing coverage/Deployment duration vs. ASV/OTU richness and number of identified species for COI and 18S per sample ###

# COI data and plots

# Format respective columns in sample_data as dates
sample_data(psCOI)<-transform(sample_data(psCOI),Deployment = as.Date(as.character(Deployment), "%Y%m%d"))
sample_data(psCOI)<-transform(sample_data(psCOI),Retrieval = as.Date(as.character(Retrieval), "%Y%m%d"))

# Calculate number of days between retrieval and deployment and add a column
sample_data(psCOI)$Deployment_Days<-sample_data(psCOI)$Retrieval-sample_data(psCOI)$Deployment

read_sums_COI<- data.frame(sample_data(psCOI)) %>% select(Observatory,Deployment_Days) # get observatories and deployment duration

sums_COI<-as.data.frame(sample_sums(psCOI)) # get read counts per sample

read_sums_COI$reads<-sums_COI[order(match(rownames(sums_COI),rownames(read_sums_COI))),]

richness_COI<-estimate_richness(psCOI,measures="Observed") # get ASV richness per sample

read_sums_COI$richness<-richness_COI[order(match(rownames(richness_COI),rownames(read_sums_COI))),]

specs_COI<- tax_glom(psCOI_species, taxrank = "Species")  # get data set agglomerated by unique species names

species_COI<-estimate_richness(specs_COI,measures="Observed") # get number of species identified in each sample

read_sums_COI$species<-species_COI[order(match(rownames(species_COI),rownames(read_sums_COI))),]

read_sums_COI$Deployment_Days<-as.numeric(read_sums_COI$Deployment_Days)

scale1=0.036 # set scale for second axis

graf_palettes$okabe_ito # get color codes of colorblind-friendly palette

read_plot_COI<-ggplot(read_sums_COI, aes(x=reads, y=richness)) +
  geom_point(aes(color="ASV/OTU richness")) +
  geom_point(aes(y = species/scale1,color="Number of species"))+
  geom_smooth(aes(reads,richness,color="ASV/OTU richness",fill="ASV/OTU richness"),method=lm)+
  geom_smooth(aes(reads,species/scale1,color="Number of species",fill="Number of species"),method=lm)+
  scale_color_manual(name="Legend",breaks=c('ASV/OTU richness', 'Number of species'),values=c('ASV/OTU richness'="#E69F00", 'Number of species'="#009E73"))+
  scale_fill_manual(name="Legend",breaks=c('ASV/OTU richness', 'Number of species'),values=c('ASV/OTU richness'="#E69F00", 'Number of species'="#009E73"))+
  scale_x_continuous(labels = scientific)+
  scale_y_continuous(limits=c(0,2500),breaks=c(0,500,1000,1500,2000,2500),labels = scales::comma,sec.axis = sec_axis(~.*scale1, name="Number of species",breaks=c(0,30,60,90)))+
  labs(x="Sequencing depth",y="ASV richness")+
  theme_bw()+
  ggtitle("C")+
  theme(legend.title = element_blank(),plot.title = element_text(size = 20,face = "bold",vjust=-.2),axis.text.y=element_text(color="black"),axis.text.x=element_text(color="black",angle=45,hjust=1),axis.text=element_text(size=10),axis.title=element_text(size=12),legend.text=element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),panel.border = element_blank())

duration_plot_COI<-ggplot(read_sums_COI, aes(x=Deployment_Days, y=richness)) +
  geom_point(aes(color="ASV/OTU richness")) +
  geom_point(aes(y = species/scale1,color="Number of species"))+
  geom_smooth(aes(Deployment_Days,richness,color="ASV/OTU richness",fill="ASV/OTU richness"),method=lm)+
  geom_smooth(aes(Deployment_Days,species/scale1,color="Number of species",fill="Number of species"),method=lm)+
  scale_color_manual(name="Legend",breaks=c('ASV/OTU richness', 'Number of species'),values=c('ASV/OTU richness'="#E69F00", 'Number of species'="#009E73"))+
  scale_fill_manual(name="Legend",breaks=c('ASV/OTU richness', 'Number of species'),values=c('ASV/OTU richness'="#E69F00", 'Number of species'="#009E73"))+
  scale_y_continuous(limits=c(0,2500),breaks=c(0,500,1000,1500,2000,2500),labels = scales::comma,sec.axis = sec_axis(~.*scale1, name="Number of species",breaks=c(0,30,60,90)))+
  scale_x_continuous(limits=c(30,NA),breaks=c(30,230,430,630))+
  labs(x="Deployment duration (days)",y="ASV richness")+
  theme_bw()+
  ggtitle("        COI\nA")+
  theme(legend.title = element_blank(),plot.title = element_text(size = 20,face = "bold",vjust=-.2),axis.text.y=element_text(color="black"),axis.text.x=element_text(color="black",angle=45,hjust=1),axis.text=element_text(size=10),axis.title=element_text(size=12),legend.text=element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),panel.border = element_blank())

# 18S data and plots

# Format respective columns in sample_data as dates
sample_data(ps18S)<-transform(sample_data(ps18S),Deployment = as.Date(as.character(Deployment), "%Y%m%d"))
sample_data(ps18S)<-transform(sample_data(ps18S),Retrieval = as.Date(as.character(Retrieval), "%Y%m%d"))

# Calculate number of days between retrieval and deployment and add a column
sample_data(ps18S)$Deployment_Days<-sample_data(ps18S)$Retrieval-sample_data(ps18S)$Deployment

read_sums_18S<- data.frame(sample_data(ps18S)) %>% select(Observatory,Deployment_Days) # get observatories and deployment duration

sums_18S<-as.data.frame(sample_sums(ps18S)) # get read counts per sample

read_sums_18S$reads<-sums_18S[order(match(rownames(sums_18S),rownames(read_sums_18S))),]

richness_18S<-estimate_richness(ps18S,measures="Observed") # get OTU richness per sample

read_sums_18S$richness<-richness_18S[order(match(rownames(richness_18S),rownames(read_sums_18S))),]

specs_18S<- tax_glom(ps18S_species, taxrank = "Species")  # get data set agglomerated by unique species names

species_18S<-estimate_richness(specs_18S,measures="Observed") # get number of species identified in each sample

read_sums_18S$species<-species_18S[order(match(rownames(species_18S),rownames(read_sums_18S))),]

read_sums_18S$Deployment_Days<-as.numeric(read_sums_18S$Deployment_Days)

scale2=(40/900) # set scale for second axis

graf_palettes$okabe_ito # get color codes of colorblind-friendly palette

read_plot_18S<-ggplot(read_sums_18S, aes(x=reads, y=richness)) +
  geom_point(aes(color="OTU richness")) +
  geom_point(aes(y = species/scale2,color="Number of species"))+
  geom_smooth(aes(reads,richness,color="OTU richness",fill="OTU richness"),method=lm)+
  geom_smooth(aes(reads,species/scale2,color="Number of species",fill="Number of species"),method=lm)+
  scale_color_manual(name="Legend",breaks=c('OTU richness', 'Number of species'),values=c('OTU richness'="#E69F00", 'Number of species'="#009E73"))+
  scale_fill_manual(name="Legend",breaks=c('OTU richness', 'Number of species'),values=c('OTU richness'="#E69F00", 'Number of species'="#009E73"))+
  scale_x_continuous(labels = scientific)+
  scale_y_continuous(limits=c(0,900),breaks=c(0,300,600,900),labels = scales::comma,sec.axis = sec_axis(~.*scale2, name="Number of species",breaks=c(0,10,20,30,40)))+
  labs(x="Sequencing depth",y="OTU richness")+
  theme_bw()+
  ggtitle("D")+
  theme(legend.title = element_blank(),plot.title = element_text(size = 20,face = "bold",vjust=-.2),axis.text.y=element_text(color="black"),axis.text.x=element_text(color="black",angle=45,hjust=1),axis.text=element_text(size=10),axis.title=element_text(size=12),legend.text=element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
       axis.line = element_line(colour = "black"),panel.border = element_blank())

duration_plot_18S<-ggplot(read_sums_18S, aes(x=Deployment_Days, y=richness)) +
  geom_point(aes(color="OTU richness")) +
  geom_point(aes(y = species/scale2,color="Number of species"))+
  geom_smooth(aes(Deployment_Days,richness,color="OTU richness",fill="OTU richness"),method=lm)+
  geom_smooth(aes(Deployment_Days,species/scale2,color="Number of species",fill="Number of species"),method=lm)+
  scale_color_manual(name="Legend",breaks=c('OTU richness', 'Number of species'),values=c('OTU richness'="#E69F00", 'Number of species'="#009E73"))+
  scale_fill_manual(name="Legend",breaks=c('OTU richness', 'Number of species'),values=c('OTU richness'="#E69F00", 'Number of species'="#009E73"))+
  scale_x_continuous(limits=c(30,NA),breaks = c(30,230,430,630))+
  scale_y_continuous(limits=c(0,900),breaks=c(0,300,600,900),labels = scales::comma,sec.axis = sec_axis(~.*scale2, name="Number of species",breaks=c(0,10,20,30,40)))+
  labs(x="Deployment duration (days)",y="OTU richness")+
  theme_bw()+
  ggtitle("        18S\nB")+
  theme(legend.title = element_blank(),plot.title = element_text(size = 20,face = "bold",vjust=-.2),axis.text.y=element_text(color="black"),axis.text.x=element_text(color="black",angle=45,hjust=1),axis.text=element_text(size=10),axis.title=element_text(size=12),legend.text=element_text(size=12),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),panel.border = element_blank())

# Combine and save plots 

read_duration_plot<-ggpubr::ggarrange(duration_plot_COI,duration_plot_18S,read_plot_COI,read_plot_18S,ncol=2,nrow=2,common.legend = T,legend="bottom",align="hv")

ggsave(read_duration_plot,file="Figures/Figure_5.png",width=6,height=7,bg="white",dpi=600)
## Outside of R, using any image processing software, move x-axis title of A and B closer to x-axes, y-axis titles of B and D closer to y-axis and plot titles of C and D closer to plots. Then save as pdf.

# Do statistics: correlation analysis and linear regression #

# COI

# Spearman Correlations
cor.test(read_sums_COI$richness,read_sums_COI$reads,method = "spearman",exact=FALSE)
cor.test(read_sums_COI$species,read_sums_COI$reads,method = "spearman",exact=FALSE)
cor.test(read_sums_COI$richness,read_sums_COI$Deployment_Days,method = "spearman",exact=FALSE)
cor.test(read_sums_COI$species,read_sums_COI$Deployment_Days,method = "spearman",exact=FALSE)

# Regression
reads_richness_COI<-lm(richness~reads,read_sums_COI)
summary(reads_richness_COI)
reads_species_COI<-lm(species~reads,read_sums_COI)
summary(reads_species_COI)

# 18S

# Spearman correlations
cor.test(read_sums_18S$richness,read_sums_18S$read,method = "spearman",exact=FALSE)
cor.test(read_sums_18S$species,read_sums_18S$reads,method = "spearman",exact=FALSE)
cor.test(read_sums_18S$richness,read_sums_18S$Deployment_Days,method = "spearman",exact=FALSE)
cor.test(read_sums_18S$species,read_sums_18S$Deployment_Days,method = "spearman",exact=FALSE)

# Regression
reads_richness_18S<-lm(richness~reads,read_sums_18S)
summary(reads_richness_18S)
reads_species_18S<-lm(species~reads,read_sums_18S)
summary(reads_species_18S)

### Number of ARMS and samples vs. ASV/OTU richness and number of identified species for COI and 18S per observatory ###

# COI data and plot

# Count samples present for each Observatory
samples1<-sample_data(psCOI)  %>% count(Observatory) # In case count() can't find the object, this is due to an issue with dplyr. Re-starting RStudio/R usually helps.
samples1<- data.frame(lapply(samples1, function(x) Reduce(c, x)))

merge_COI<-merge_samples(psCOI,"Observatory") # merge Samples for each Observatory

sums_merge_COI<-as.data.frame(sample_names(merge_COI))
colnames(sums_merge_COI)<-"Observatory"

sums_merge_COI$samples<-samples1$n # Add number of samples

sums_merge_COI$arms<-c(2,2,5,4,2,3,1,6,3,5,5,8,4,2,4) # Manually add number of ARMS units per Observatory

read_sums_merge_COI<-as.data.frame(sample_sums(merge_COI)) # get read counts per Observatory

sums_merge_COI$reads<-read_sums_merge_COI[order(match(rownames(read_sums_merge_COI),rownames(sums_merge_COI))),]

richness_merge_COI<-estimate_richness(merge_COI,measures="Observed") # get ASV richness per Observatory

sums_merge_COI$richness<-richness_merge_COI[order(match(rownames(richness_merge_COI),rownames(sums_merge_COI))),]

merge_COI_species<-subset_taxa(merge_COI,!is.na(Species))

specs_merge_COI<- tax_glom(merge_COI_species, taxrank = "Species")  # get data set agglomerated by unique species names

species_merge_COI<-estimate_richness(specs_merge_COI,measures="Observed") # get number of species identified in each observatory

sums_merge_COI$species<-species_merge_COI[order(match(rownames(species_merge_COI),rownames(sums_merge_COI))),]

scale1 <- 25 # set scale for plot

sums_merge_COI<-sums_merge_COI %>% mutate(species = species*scale1) 

# Make labels for barplots
sums_merge_COI$label_richness<-ifelse(sums_merge_COI$richness > sums_merge_COI$species,paste0("n=",as.character(sums_merge_COI$arms),"\n(",as.character(sums_merge_COI$samples),")"),"")
sums_merge_COI$label_species<-ifelse(sums_merge_COI$richness < sums_merge_COI$species,paste0("n=",as.character(sums_merge_COI$arms),"\n(",as.character(sums_merge_COI$samples),")"),"")

sums_merge_COI<-sums_merge_COI %>% pivot_longer(c(richness,species))
labels1<-sums_merge_COI %>% select(-c(name,value)) %>% pivot_longer(c(label_richness,label_species)) # Label dataframe for barplot
labels1<-labels1 %>% distinct()

sums_merge_COI<-sums_merge_COI[order(sums_merge_COI$arms),] # Order by number of ARMS units for plot
labels1<-labels1[order(labels1$arms),] # Order by number of ARMs units for plot

reads_line_COI <- sums_merge_COI %>% select(c(Observatory,reads)) 

reads_line_COI<-reads_line_COI %>% mutate(reads = reads/50) 

observatory_plot_COI<-ggplot(sums_merge_COI, aes(x=Observatory)) +
  geom_bar( aes(y = value, fill = name, group = name),stat="identity", position=position_dodge(),color="black") +
  geom_line(data=reads_line_COI, aes(x=Observatory, y=reads, group=1))+
  scale_fill_manual(values = c("#E69F00", "#009E73"),labels=c('ASV/OTU richness', 'Number of species')) +
  scale_y_continuous(limits=c(0,7900),breaks=c(0,2500,5000,7500),name = expression(atop("ASV richness", "Sequencing depth [x50]")),labels = scales::comma,sec.axis = sec_axis(~./scale1, name="Number of species",labels = scales::comma,breaks=c(0,100,200,300)))+
  scale_x_discrete(expand=expansion(mult = c(0, 0),add = c(1.5,.5)),limits=sums_merge_COI$Observatory)+
  theme_bw()+
  ggtitle("A")+
  theme(legend.key.size = unit(.75,"line"),legend.text=element_text(size=10),axis.title=element_text(size=11),axis.ticks.x = element_blank(),axis.title.x=element_blank(),legend.title = element_blank(),plot.title = element_text(size = 20,face = "bold",vjust=-.2),axis.text.y=element_text(color="black"),axis.text.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),panel.border = element_blank(),panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),legend.background = element_rect(color = NA,fill='transparent'))+
  geom_text(aes(label=labels1$value,y=value,group = Observatory), vjust = -0.5,position = position_dodge(width = 0.9),size=3) 

# 18S data and plot

# Count samples present for each Observatory
samples2<-sample_data(ps18S)  %>% count(Observatory) # In case count() can't find the object, this is due to an issue with dplyr. Re-starting RStudio/R usually helps.
samples2<- data.frame(lapply(samples2, function(x) Reduce(c, x)))

merge_18S<-merge_samples(ps18S,"Observatory") # merge Samples for each Observatory

sums_merge_18S<-as.data.frame(sample_names(merge_18S))
colnames(sums_merge_18S)<-"Observatory"

sums_merge_18S$samples<-samples2$n # Add number of samples

sums_merge_18S$arms<-c(2,2,5,4,2,3,1,6,3,5,5,8,4,2,4) # Manually add number of ARMS units per Observatory

read_sums_merge_18S<-as.data.frame(sample_sums(merge_18S)) # get read counts per Observatory

sums_merge_18S$reads<-read_sums_merge_18S[order(match(rownames(read_sums_merge_18S),rownames(sums_merge_18S))),]

richness_merge_18S<-estimate_richness(merge_18S,measures="Observed") # get ASV richness per Observatory

sums_merge_18S$richness<-richness_merge_18S[order(match(rownames(richness_merge_18S),rownames(sums_merge_18S))),]

merge_18S_species<-subset_taxa(merge_18S,!is.na(Species))

specs_merge_18S<- tax_glom(merge_18S_species, taxrank = "Species")  # get data set agglomerated by unique species names

species_merge_18S<-estimate_richness(specs_merge_18S,measures="Observed") # get number of species identified in each observatory

sums_merge_18S$species<-species_merge_18S[order(match(rownames(species_merge_18S),rownames(sums_merge_18S))),]

scale2 <- 25 # set scale for plot

sums_merge_18S<-sums_merge_18S %>% mutate(species = species*scale2) 

# Make labels for barplots
sums_merge_18S$label_richness<-ifelse(sums_merge_18S$richness > sums_merge_18S$species,paste0("n=",as.character(sums_merge_18S$arms),"\n(",as.character(sums_merge_18S$samples),")"),"")
sums_merge_18S$label_species<-ifelse(sums_merge_18S$richness < sums_merge_18S$species,paste0("n=",as.character(sums_merge_18S$arms),"\n(",as.character(sums_merge_18S$samples),")"),"")

sums_merge_18S<-sums_merge_18S %>% pivot_longer(c(richness,species))
labels2<-sums_merge_18S %>% select(-c(name,value)) %>% pivot_longer(c(label_richness,label_species)) # Label dataframe for barplot
labels2<-labels2 %>% distinct()

sums_merge_18S<-sums_merge_18S[order(sums_merge_18S$arms),] # Order by number of ARMS units for plot
labels2<-labels2[order(labels2$arms),] # Order by number of ARMS units for plot

reads_line_18S <- sums_merge_18S %>% select(c(Observatory,reads)) 

reads_line_18S<-reads_line_18S %>% mutate(reads = reads/400) 

observatory_plot_18S<-ggplot(sums_merge_18S, aes(x=Observatory)) +
  geom_bar( aes(y = value, fill = name, group = name),stat="identity", position=position_dodge(),color="black") +
  geom_line(data=reads_line_18S, aes(x=Observatory, y=reads, group=1))+
  scale_fill_manual(values = c("#E69F00", "#009E73"),labels=c('ASV/OTU richness', 'Number of species')) +
  scale_y_continuous(limits=c(0,2700),breaks=c(0,500,1000,1500,2000,2500),name = expression(atop("OTU richness", "Sequencing depth [x400]")),labels = scales::comma,sec.axis = sec_axis(~./scale1, name="Number of species",labels = scales::comma,breaks=c(0,50,100,150,200,250)))+
  scale_x_discrete(expand=expansion(mult = c(0, 0),add = c(1.5,.5)),limits=sums_merge_18S$Observatory)+
  theme_bw()+
  ggtitle("B")+
  theme(legend.key.size = unit(.75,"line"),legend.text=element_text(size=10),axis.title=element_text(size=11),axis.title.x=element_blank(),legend.title = element_blank(),plot.title = element_text(size = 20,face = "bold",vjust=-.2),axis.text.y=element_text(color="black"),axis.text.x=element_text(color="black",angle=45,hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),panel.border = element_blank(),panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA),legend.background = element_rect(color = NA,fill='transparent'))+
  geom_text(aes(label=labels2$value,y=value,group = Observatory), vjust = -0.5,position = position_dodge(width = 0.9),size=3) 

# Combine plots

observatory_plot<-ggpubr::ggarrange(observatory_plot_COI,observatory_plot_18S,ncol=1,nrow=2,common.legend = T,legend="bottom",align="hv")
ggsave(observatory_plot,file="Figures/Figure_3.png",width=5,height=6,bg="transparent",dpi=600)
## Outside of R, using any image processing software, adjust space between plots, add arrow and accompanying text. ##

# Make scatter plot for diversity vs. sampling effort per observatory

sums_merge_COI_wide<-sums_merge_COI %>% pivot_wider(names_from = name)
sums_merge_COI_wide<-sums_merge_COI_wide %>% mutate(species = species/scale1) 

sums_merge_18S_wide<-sums_merge_18S %>% pivot_wider(names_from = name)
sums_merge_18S_wide<-sums_merge_18S_wide %>% mutate(species = species/scale2)

scale3=(300/6000)

arms_plot_COI<-ggplot(sums_merge_COI_wide, aes(x=arms, y=richness)) +
  geom_point(aes(color="ASV/OTU richness")) +
  geom_point(aes(y = species/scale3,color="Number of species"))+
  geom_smooth(aes(arms,richness,color="ASV/OTU richness",fill="ASV/OTU richness"),method=lm)+
  geom_smooth(aes(arms,species/scale3,color="Number of species",fill="Number of species"),method=lm)+
  scale_color_manual(name="Legend",breaks=c('ASV/OTU richness', 'Number of species'),values=c('ASV/OTU richness'="#E69F00", 'Number of species'="#009E73"))+
  scale_fill_manual(name="Legend",breaks=c('ASV/OTU richness', 'Number of species'),values=c('ASV/OTU richness'="#E69F00", 'Number of species'="#009E73"))+
  scale_y_continuous(limits=c(0,6700),breaks=c(0,2000,4000,6000),labels = scales::comma,sec.axis = sec_axis(~.*scale3, name="Number of species",breaks=c(0,100,200,300)))+
  scale_x_continuous(limits = c(1,8),breaks=1:9)+
  labs(x="Number of ARMS units",y="ASV richness")+
  theme_bw()+
  ggtitle("               COI\nA")+
  theme(legend.title = element_blank(),plot.title = element_text(size = 20,face = "bold",vjust=-.2),axis.text.y=element_text(color="black"),axis.text.x=element_text(color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),panel.border = element_blank())

samples_plot_COI<-ggplot(sums_merge_COI_wide, aes(x=samples, y=richness)) +
  geom_point(aes(color="ASV/OTU richness")) +
  geom_point(aes(y = species/scale3,color="Number of species"))+
  geom_smooth(aes(samples,richness,color="ASV/OTU richness",fill="ASV/OTU richness"),method=lm)+
  geom_smooth(aes(samples,species/scale3,color="Number of species",fill="Number of species"),method=lm)+
  scale_color_manual(name="Legend",breaks=c('ASV/OTU richness', 'Number of species'),values=c('ASV/OTU richness'="#E69F00", 'Number of species'="#009E73"))+
  scale_fill_manual(name="Legend",breaks=c('ASV/OTU richness', 'Number of species'),values=c('ASV/OTU richness'="#E69F00", 'Number of species'="#009E73"))+
  scale_y_continuous(limits=c(0,6700),breaks=c(0,2000,4000,6000),labels = scales::comma,sec.axis = sec_axis(~.*scale3, name="Number of species",breaks=c(0,100,200,300)))+
  scale_x_continuous(limits = c(3,24),breaks=seq(3,24,3))+
  labs(x="Number of samples",y="ASV richness")+
  theme_bw()+
  ggtitle("\nC")+
  theme(legend.title = element_blank(),plot.title = element_text(size = 20,face = "bold",vjust=-.2),axis.text.y=element_text(color="black"),axis.text.x=element_text(color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),panel.border = element_blank())

scale4=(75/2500)

arms_plot_18S<-ggplot(sums_merge_18S_wide, aes(x=arms, y=richness)) +
  geom_point(aes(color="ASV/OTU richness")) +
  geom_point(aes(y = species/scale4,color="Number of species"))+
  geom_smooth(aes(arms,richness,color="ASV/OTU richness",fill="ASV/OTU richness"),method=lm)+
  geom_smooth(aes(arms,species/scale4,color="Number of species",fill="Number of species"),method=lm)+
  scale_color_manual(name="Legend",breaks=c('ASV/OTU richness', 'Number of species'),values=c('ASV/OTU richness'="#E69F00", 'Number of species'="#009E73"))+
  scale_fill_manual(name="Legend",breaks=c('ASV/OTU richness', 'Number of species'),values=c('ASV/OTU richness'="#E69F00", 'Number of species'="#009E73"))+
  scale_y_continuous(limits=c(0,2500),breaks=c(0,500,1000,1500,2000,2500),labels = scales::comma,sec.axis = sec_axis(~.*scale4, name="Number of species",breaks=c(0,25,50,75)))+
  scale_x_continuous(limits = c(1,8),breaks=1:8)+
  labs(x="Number of ARMS units",y="OTU richness")+
  theme_bw()+
  ggtitle("               18S\nB")+
  theme(legend.title = element_blank(),plot.title = element_text(size = 20,face = "bold",vjust=-.2),axis.text.y=element_text(color="black"),axis.text.x=element_text(color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),panel.border = element_blank())

samples_plot_18S<-ggplot(sums_merge_18S_wide, aes(x=samples, y=richness)) +
  geom_point(aes(color="ASV/OTU richness")) +
  geom_point(aes(y = species/scale4,color="Number of species"))+
  geom_smooth(aes(samples,richness,color="ASV/OTU richness",fill="ASV/OTU richness"),method=lm)+
  geom_smooth(aes(samples,species/scale4,color="Number of species",fill="Number of species"),method=lm)+
  scale_color_manual(name="Legend",breaks=c('ASV/OTU richness', 'Number of species'),values=c('ASV/OTU richness'="#E69F00", 'Number of species'="#009E73"))+
  scale_fill_manual(name="Legend",breaks=c('ASV/OTU richness', 'Number of species'),values=c('ASV/OTU richness'="#E69F00", 'Number of species'="#009E73"))+
  scale_y_continuous(limits=c(0,2500),breaks=c(0,500,1000,1500,2000,2500),labels = scales::comma,sec.axis = sec_axis(~.*scale4, name="Number of species",breaks=c(0,25,50,75)))+
  scale_x_continuous(limits = c(3,24),breaks=seq(3,24,3))+
  labs(x="Number of samples",y="OTU richness")+
  theme_bw()+
  ggtitle("\nD")+
  theme(legend.title = element_blank(),plot.title = element_text(size = 20,face = "bold",vjust=-.2),axis.text.y=element_text(color="black"),axis.text.x=element_text(color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),panel.border = element_blank())

# Combine and save plots 

effort_plot<-ggpubr::ggarrange(arms_plot_COI,arms_plot_18S,samples_plot_COI,samples_plot_18S,ncol=2,nrow=2,common.legend = T,legend="bottom",align="hv")

ggsave(effort_plot,file="effort_plot.png",width=8,height=8,bg="white")

# Do statistics: correlation analysis and linear regression #

# COI

# Pearson Correlations

cor.test(sums_merge_COI_wide$species,sums_merge_COI_wide$arms,method = "spearman",exact=FALSE)
cor.test(sums_merge_COI_wide$richness,sums_merge_COI_wide$arms,method = "spearman",exact=FALSE)
cor.test(sums_merge_COI_wide$species,sums_merge_COI_wide$samples,method = "spearman",exact=FALSE)
cor.test(sums_merge_COI_wide$richness,sums_merge_COI_wide$samples,method = "spearman",exact=FALSE)

# Regression
effort1_species_COI<-lm(species~arms,sums_merge_COI_wide)
summary(effort1_species_COI)
effort2_species_COI<-lm(species~samples,sums_merge_COI_wide)
summary(effort2_species_COI)

# 18S

# Pearson correlations
cor.test(sums_merge_18S_wide$species,sums_merge_18S_wide$arms,method = "spearman",exact=FALSE)
cor.test(sums_merge_18S_wide$richness,sums_merge_18S_wide$arms,method = "spearman",exact=FALSE)
cor.test(sums_merge_18S_wide$species,sums_merge_18S_wide$samples,method = "spearman",exact=FALSE)
cor.test(sums_merge_18S_wide$richness,sums_merge_18S_wide$samples,method = "spearman",exact=FALSE)

# Regression
effort1_species_18S<-lm(species~arms,sums_merge_18S_wide)
summary(effort1_species_18S)
effort2_species_18S<-lm(species~samples,sums_merge_18S_wide)
summary(effort2_species_18S)
effort1_richness_18S<-lm(richness~arms,sums_merge_18S_wide)
summary(effort1_richness_18S)
effort2_richness_18S<-lm(richness~samples,sums_merge_18S_wide)
summary(effort2_richness_18S)

# Save info on number of ASVs/OTUs and identified species per observatory ##

write.table(sums_merge_COI_wide[,-(5:6)],"observatory_ASV_species_richness_COI.txt",sep="\t",row.names = F)
write.table(sums_merge_18S_wide[,-(5:6)],"observatory_OTU_species_richness_18S.txt",sep="\t",row.names = F)

## Make frequency plots for ASVs/OTUs and number of identified species per observatory ##

# COI ASV and 18S OTU frequency tables

asv_freq<-as.data.frame(rowSums(as.data.frame(t(otu_table(merge_COI)))!=0))
colnames(asv_freq)<-"freq"
asv_freq<-asv_freq%>% group_by(freq) %>% tally()

otu_freq<-as.data.frame(rowSums(as.data.frame(t(otu_table(merge_18S)))!=0))
colnames(otu_freq)<-"freq"
otu_freq<-otu_freq%>% group_by(freq) %>% tally()

# COI and 18S identified species frequency tables

coi_spec_freq<-as.data.frame(rowSums(as.data.frame(t(otu_table(specs_merge_COI)))!=0))
colnames(coi_spec_freq)<-"freq"
coi_spec_freq<-coi_spec_freq%>% group_by(freq) %>% tally()

specs_18S_freq<-as.data.frame(rowSums(as.data.frame(t(otu_table(specs_merge_18S)))!=0))
colnames(specs_18S_freq)<-"freq"
specs_18S_freq<-specs_18S_freq%>% group_by(freq) %>% tally()

asv_freq_plot<-ggplot(asv_freq, aes(x=freq,y=n)) +
  geom_bar(stat="identity", fill="#E69F00",color="black") +
  theme_bw()+
  scale_y_continuous(labels = scales::comma,trans="log2",breaks=c(1,10,100,1000,7000,40000))+
  scale_x_continuous(breaks=c(1,4,7,10))+
  ggtitle("           COI\nA")+
  theme(plot.margin=unit(c(0,0,0,0), "cm"),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.title = element_blank(),plot.title = element_text(size = 16,face = "bold",vjust=-.2),axis.text.y=element_text(color="black"),axis.text.x=element_text(color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),panel.border = element_blank(),panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA))

otu_freq_plot<-ggplot(otu_freq, aes(x=freq,y=n))+
  geom_bar(stat="identity", fill="#E69F00",color="black") +
  theme_bw()+
  scale_y_continuous(labels = scales::comma,trans="log2",breaks=c(1,10,100,1000,6000))+
  scale_x_continuous(breaks=c(1,4,7,10,13))+
  ggtitle("           18S\nB")+
  theme(plot.margin=unit(c(0,0,0,-3), "cm"),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.title = element_blank(),plot.title = element_text(size = 16,face = "bold",vjust=-.2),axis.text.y=element_text(color="black"),axis.text.x=element_text(color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),panel.border = element_blank(),panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA))

coi_spec_freq_plot<-ggplot(coi_spec_freq, aes(x=freq,y=n)) +
  geom_bar(stat="identity", fill="#009E73",color="black") +
  theme_bw()+
  scale_y_continuous(labels = scales::comma)+
  scale_x_continuous(breaks=c(1,4,7,10,13))+
  ggtitle("\nC")+
  theme(plot.margin=unit(c(0,0,0,0), "cm"),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.title = element_blank(),plot.title = element_text(size = 16,face = "bold",vjust=-.2),axis.text.y=element_text(color="black"),axis.text.x=element_text(color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),panel.border = element_blank(),panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA))

specs_18S_freq_plot<-ggplot(specs_18S_freq, aes(x=freq,y=n)) +
  geom_bar(stat="identity", fill="#009E73",color="black") +
  theme_bw()+
  scale_y_continuous(labels = scales::comma,limits=c(0,75),breaks=c(0,25,50,75))+
  scale_x_continuous(breaks=c(1,4,7,10))+
  ggtitle("\nD")+
  theme(plot.margin=unit(c(0,0,0,-3), "cm"),axis.title.x=element_blank(),axis.title.y=element_blank(),legend.title = element_blank(),plot.title = element_text(size = 16,face = "bold",vjust=-.2),axis.text.y=element_text(color="black"),axis.text.x=element_text(color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),panel.border = element_blank(),panel.background = element_rect(fill='transparent'),plot.background = element_rect(fill='transparent', color=NA))

# Combine and save plots 

freq_plot<-ggarrange(asv_freq_plot,otu_freq_plot,coi_spec_freq_plot,specs_18S_freq_plot,ncol=2,nrow=2,align="hv")
freq_plot<-annotate_figure(freq_plot,bottom = textGrob("Frequency",gp = gpar(fontsize = 12)),left = textGrob("Count", gp = gpar(fontsize = 12),rot=90))

ggsave(freq_plot,file="Figures/Figure_4.png",width=5,height=4,bg="transparent",dpi=600)
# Use legend of Figure_3.png, add it to this plot outside of R using any image processing software. Also, move plots closer together.

# Calculate some values for frequency distributions

38090/sum(asv_freq$n) # Share of COI ASVs occurring at one observatory only

1593/sum(asv_freq$n) # Share of COI ASVs occurring at two observatories only

5794/sum(otu_freq$n) # Share of 18S OTUs occurring at one observatory only

1565/sum(otu_freq$n) # Share of 18S OTUs occurring at two observatories only

393/sum(coi_spec_freq$n) # Share of Species occurring at one observatory only for COI

160/sum(coi_spec_freq$n) # Share of Species occurring at two observatories only for COI

13/sum(coi_spec_freq$n) # Share of Species occurring at 10 or more observatories for COI

68/sum(specs_18S_freq$n) # Share of Species occurring at one observatory only for 18S

17/sum(specs_18S_freq$n) # Share of Species occurring at two observatories only for 18S

2/sum(specs_18S_freq$n) # Share of Species occurring at 10 or more observatories for 18S

## Compare diversity among samples with varying anthropogenic influence types ##

# Remove samples with a read number of < 5,000

psCOI_reads <- prune_samples(sample_sums(psCOI) > 5000, psCOI)
ps18S_reads <- prune_samples(sample_sums(ps18S) > 5000, ps18S)

# rarefy samples to 5,000 reads ( a large number of ASVs/OTUs will be removed, as the previous step left many sequences with zero reads in the data sets)

psCOI_rarefied <- rarefy_even_depth(psCOI_reads, rngseed=1, sample.size=5000, replace=F)
ps18S_rarefied <- rarefy_even_depth(ps18S_reads, rngseed=1, sample.size=5000, replace=F)

# Make phyloseq objects for ASVs/OTUs assigned to species level only and agglomerate to unique species

psCOI_rarefied_spec<-subset_taxa(psCOI_rarefied,!is.na(Species))
psCOI_rarefied_spec<-tax_glom(psCOI_rarefied_spec,"Species")

ps18S_rarefied_spec<-subset_taxa(ps18S_rarefied,!is.na(Species))
ps18S_rarefied_spec<-tax_glom(ps18S_rarefied_spec,"Species")

# Get ASV/OTU and species richness for samples and combine tables

COI_rarefied_richness<-estimate_richness(psCOI_rarefied,measure="Observed")
colnames(COI_rarefied_richness)<-"seq_rich"

rarefied_richness_18S<-estimate_richness(ps18S_rarefied,measure="Observed")
colnames(rarefied_richness_18S)<-"seq_rich"

COI_rarefied_spec_richness<-estimate_richness(psCOI_rarefied_spec,measure="Observed")
colnames(COI_rarefied_spec_richness)<-"spec_rich"

rarefied_spec_richness_18S<-estimate_richness(ps18S_rarefied_spec,measure="Observed")
colnames(rarefied_spec_richness_18S)<-"spec_rich"

COI_rich<-cbind(COI_rarefied_richness,COI_rarefied_spec_richness)
rich_18S<-cbind(rarefied_richness_18S,rarefied_spec_richness_18S)

# Add influence info

Anthropogenic_influence<-cbind(rownames(samples),samples$Anthropogenic_influence)
colnames(Anthropogenic_influence)<-c("Sample","Anthropogenic_influence")

COI_rich$Sample<-rownames(COI_rich)
COI_rich<-merge(COI_rich,Anthropogenic_influence,by="Sample")

rich_18S$Sample<-rownames(rich_18S)
rich_18S<-merge(rich_18S,Anthropogenic_influence,by="Sample")

# Make Industrial and Semi-industrial same group, as there are only a handful of Industrial samples 
COI_rich$Anthropogenic_influence<-ifelse(COI_rich$Anthropogenic_influence=="industrial" | COI_rich$Anthropogenic_influence=="semi-industrial","industrial/semi-industrial",COI_rich$Anthropogenic_influence)

rich_18S$Anthropogenic_influence<-ifelse(rich_18S$Anthropogenic_influence=="industrial" | rich_18S$Anthropogenic_influence=="semi-industrial","industrial/semi-industrial",rich_18S$Anthropogenic_influence)

# Mean & SD

as.data.frame(COI_rich %>% # as.data.frame necessary, cause RStudio does not display decimal digits properly for tibbles
  group_by(Anthropogenic_influence) %>%
  summarise_at(vars(seq_rich), list(name = mean)))

as.data.frame(COI_rich %>% 
                group_by(Anthropogenic_influence) %>%
                summarise_at(vars(seq_rich), list(name = sd)))

as.data.frame(COI_rich %>% 
                group_by(Anthropogenic_influence) %>%
                summarise_at(vars(spec_rich), list(name = mean)))

as.data.frame(COI_rich %>% 
                group_by(Anthropogenic_influence) %>%
                summarise_at(vars(spec_rich), list(name = sd)))

as.data.frame(rich_18S %>% 
                group_by(Anthropogenic_influence) %>%
                summarise_at(vars(seq_rich), list(name = mean)))

as.data.frame(rich_18S %>% 
                group_by(Anthropogenic_influence) %>%
                summarise_at(vars(seq_rich), list(name = sd)))

as.data.frame(rich_18S %>% 
                group_by(Anthropogenic_influence) %>%
                summarise_at(vars(spec_rich), list(name = mean)))

as.data.frame(rich_18S %>% 
                group_by(Anthropogenic_influence) %>%
                summarise_at(vars(spec_rich), list(name = sd)))

# Check for normality and statistically compare ASV/OTU and species richness between Anthropogenic_influences

shapiro.test(COI_rich$seq_rich) # Significant
shapiro.test(log1p(COI_rich$seq_rich)) # Not significant

shapiro.test(COI_rich$spec_rich) # Not Significant

COI_rich$seq_rich_log<-log1p(COI_rich$seq_rich) # log(1+x)-transform
aov_COI_seq_rich<-aov(seq_rich_log ~ Anthropogenic_influence,COI_rich)
summary(aov_COI_seq_rich) # not significant

aov_COI_spec_rich<-aov(spec_rich ~ Anthropogenic_influence,COI_rich)
summary(aov_COI_spec_rich) # significant, do post-hoc test
TukeyHSD(aov_COI_spec_rich)

shapiro.test(rich_18S$seq_rich) # Not significant

shapiro.test(rich_18S$spec_rich) # Significant
shapiro.test(log1p(rich_18S$spec_rich)) # Significant

aov_18S_seq_rich<-aov(seq_rich ~ Anthropogenic_influence,rich_18S)
summary(aov_18S_seq_rich) # not significant

kruskal.test(spec_rich~Anthropogenic_influence,rich_18S) 

# Plots

COI_rich<-COI_rich %>% select(-seq_rich_log)
COI_rich$spec_rich<-COI_rich$spec_rich*10
COI_rich<-COI_rich %>% pivot_longer(cols=-c(Sample,Anthropogenic_influence),names_to = "type", values_to = "richness")
COI_rich$Anthropogenic_influence<-gsub("industrial/semi-industrial","(semi)industrial",COI_rich$Anthropogenic_influence)
COI_rich$Anthropogenic_influence<-gsub(".*LHI.*","LHI",COI_rich$Anthropogenic_influence)
COI_rich$type<-as.factor(COI_rich$type)

coi_plot<-ggplot(COI_rich, aes(x = Anthropogenic_influence)) +
  geom_boxplot(aes(y = richness, fill = type),outlier.size=.5)+ 
  theme_classic()+
  scale_y_continuous(name="ASV richness",sec.axis = sec_axis(~.*.1, name="Number of species"),limits=c(0,820))+
  scale_fill_manual(values=c(seq_rich="#E69F00", spec_rich="#009E73"),labels=c("ASV/OTU richness","Number of species"))+
  theme(legend.title = element_blank(),plot.title = element_text(size = 16,face = "bold",vjust=-.2),axis.title.x = element_blank(),axis.text.y=element_text(color="black"),panel.background = element_rect(fill = "transparent",colour = NA_character_),plot.background = element_rect(fill = "transparent",colour = NA_character_),legend.position="none",axis.text.x=element_text(angle=45,hjust=1,color="black"))+
  ggtitle("A")+
  labs(x="Anthropogenic_influence type",y="ASV richness")

rich_18S$spec_rich<-rich_18S$spec_rich*20
rich_18S<-rich_18S %>% pivot_longer(cols=-c(Sample,Anthropogenic_influence),names_to = "type", values_to = "richness")
rich_18S$Anthropogenic_influence<-gsub("industrial/semi-industrial","(semi)industrial",rich_18S$Anthropogenic_influence)
rich_18S$Anthropogenic_influence<-gsub(".*LHI.*","LHI",rich_18S$Anthropogenic_influence)
rich_18S$type<-as.factor(rich_18S$type)

plot_18S<-ggplot(rich_18S, aes(x = Anthropogenic_influence)) +
  geom_boxplot(aes(y = richness, fill = type),outlier.size=.5)+ 
  theme_classic()+
  scale_y_continuous(name="OTU richness",sec.axis = sec_axis(~.*.05, name="Number of species"))+
  scale_fill_manual(values=c(seq_rich="#E69F00", spec_rich="#009E73"),labels=c("ASV/OTU richness","Number of species"))+
  theme(legend.title = element_blank(),plot.title = element_text(size = 16,face = "bold",vjust=-.2),axis.title.x = element_blank(),axis.text.y=element_text(color="black"),panel.background = element_rect(fill = "transparent",colour = NA_character_),plot.background = element_rect(fill = "transparent",colour = NA_character_),legend.position="none",axis.text.x=element_text(angle=45,hjust=1,color="black"))+
  ggtitle("B")+
  labs(x="Anthropogenic_influence type",y="OTU richness")

richness_plot<-ggpubr::ggarrange(coi_plot,plot_18S,ncol=2,nrow=1,common.legend = T,legend="bottom",align="hv")

ggsave(richness_plot,file="Figures/Figure_7.png",width=5,height=3.5,bg="white")  
# Draw significance denotions for COI number of species using any image processing software outside of R (geom_signif did not work here with grouped boxplot)

