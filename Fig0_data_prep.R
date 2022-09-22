##############################################################################################################################################
#+ 1) Purposes: conduct correlation analysis among the main exposure and outcome variables, 
#+              Main exposure: omega3 - long/short term, ala epa dha dpa
#+ 2) Study design:  Cross-sectional design
#+ 3) Endpoints: primary: 
# 4) Exposures: primary: cumulative average of long-term physical activity level
# 5) Covariates: age, diet quality, total energy intake, smoking status, antibiotic use, probiotic use, stool type
# 6) Follow-up: 2012 (Men's Lifestyle Validation Study)  
#############################################################################################################################################

rm(list=ls())
getwd()
setwd("./noGit")
options(max.print=1000000)

library(viridis)
library(ggplot2)
library(grid)
library(tidyverse)
library(cowplot)
library(data.table)
library(readr)
library(haven)
library(plyr)
library(dplyr)
library(vegan)
library(scales)
library(RColorBrewer)
library(grid)
library(pheatmap)
library(lme4)
library(ggplot2)
library(cowplot)
library(ggsignif)
library(car)
library(Maaslin2)
library(gridExtra)
library(knitr)
library(readxl)
library(ggpubr)
library(tidyverse)
library(vegan)
library(RColorBrewer)
library(grid)
library(ggpubr)
library(viridis)
library("Hmisc")
library(reshape2)
library(expss)
source("Rstart.R")


###### DATA PREPARATION ########

#### Import Metadata by individual #####
meta_id <- read.table("data/mlvs_exposure_for_jorick.txt", header=TRUE, sep='\t', check.names=TRUE, quote ="")
idkey <- read.csv('idkey.csv', header=TRUE, sep = ",", stringsAsFactors = FALSE)
rownames(idkey) <- idkey$id
meta_id<- left_join(idkey, meta_id, by="id")

#### Import Plasma data by individual#####
plasma <- read.table("data/plasma_fatty_acids_sjbae.txt",
                       header=TRUE, sep='\t', check.names=TRUE, quote ="")
plasma$ID1 <- as.integer(substr(plasma$ID, 1, 6)) 
#get first 6 number of ID1 which stands for id in idkey file

#### Merge and get individual metadata 
meta_pl_id<- inner_join(meta_id, plasma, by ="ID1")

#### set weeks #####
colnames(meta_pl_id) <- gsub("w1avg","w1", colnames(meta_pl_id)) 
colnames(meta_pl_id) <- gsub("w2avg","w2", colnames(meta_pl_id))
colnames(meta_pl_id) <- gsub("plasma1","w1", colnames(meta_pl_id))
colnames(meta_pl_id) <- gsub("plasma2","w2", colnames(meta_pl_id))

meta_pl_id$ala_w1 <- meta_pl_id$a_ala_fs_dr_w1
meta_pl_id$epa_w1 <- meta_pl_id$a_f205_fs_dr_w1
meta_pl_id$dha_w1 <- meta_pl_id$a_f226_fs_dr_w1
meta_pl_id$dpa_w1 <- meta_pl_id$a_p22_5_fs_dr_w1
meta_pl_id$omega3_noala_w1 <- meta_pl_id$a_omega3_noala_fs_dr_w1
meta_pl_id$omega3_w1 <- meta_pl_id$a_omega3_fs_dr_w1

meta_pl_id$ala_w2 <- meta_pl_id$a_ala_fs_dr_w2
meta_pl_id$epa_w2 <- meta_pl_id$a_f205_fs_dr_w2
meta_pl_id$dha_w2 <- meta_pl_id$a_f226_fs_dr_w2
meta_pl_id$dpa_w2 <- meta_pl_id$a_p22_5_fs_dr_w2
meta_pl_id$omega3_noala_w2 <- meta_pl_id$a_omega3_noala_fs_dr_w2
meta_pl_id$omega3_w2 <- meta_pl_id$a_omega3_fs_dr_w2

##### save the metadata of exposure and plasma outcome ####

write.table(meta_pl_id, file = "meta_pl_id.tsv", sep = "\t", col.names = TRUE, row.names = TRUE)


#### TAXONOMY DATA ####
# read in taxonomy data
# only keep species-level features

all_species <- read_tsv("data/metaphlan_taxonomic_profiles.tsv")
names(all_species)[names(all_species) == '# taxonomy'] <- 'Sample'

all_species <-all_species %>%
  separate(Sample, c("kingdom","phylum","class" ,"order","family","genus" ,"species" ), 
           sep = '\\|', remove = TRUE)

all_species1 <- subset(all_species,!is.na(species))
all_species1 <- as.data.frame(all_species1) #weird formatting change to maintain rows
rownames(all_species1)<-all_species1$species

###### TAX and OTU table
TAX <- as.matrix(all_species1[,1:7])
OTU <- all_species1[,8:ncol(all_species1)]


all_species1<-all_species1[,-c(1:7)]
all_species1 <- as.data.frame(all_species1)
dim(all_species1)
# sunjeong: [1] 543 929

all_sample <- colnames(all_species1)
#View(all_sample)

# remove guy with colectomy 
all_species1 <- all_species1[ , !colnames(all_species1) %in% c(grep("15074981", colnames(all_species1), value = T))]


# generate list of species after abundance/prevalence filtering
# since this file is scaled to 100, i am keeping bugs with >10% prevalence and 0.01 abundance (i.e. .0001*100 when relab is scaled 0-1 and not 0-100)
dim(all_species1[apply(all_species1, 1, function(bug) sum(bug >= 0.01) >= 0.1*ncol(all_species1)), ])
# [1] 139 925 = the 139 species after QC included in initial starr manuscripts
# sunjeong: [1] 156 925
all_species_filt <- all_species1[ apply(all_species1, 1, function(bug) sum(bug >= 0.01) >= 0.1*ncol(all_species1)), ]


# important! bug tables from channing need to be converted to 0-1 (not 0-11), which will wreck our downstream stats
all_species_filt   <- all_species_filt  / 100 

all_species_filt <- t(all_species_filt)
all_species_filt <- as.data.frame(all_species_filt)

all_species_filt$ID1 <- substr(rownames(all_species_filt), start = 1, stop = 6)
all_species_filt$stn <- substr(rownames(all_species_filt), start = 10, stop = 13)
all_species_filt$stn <- str_replace_all(all_species_filt$stn,  c("SF05" = "s1","SF06" = "s2","SF07" = "s3","SF08" = "s4"))
rownames(all_species_filt) <- paste0(all_species_filt$ID1,"_", all_species_filt$stn)

species_filt <- all_species_filt[,!colnames(all_species_filt) %in% c('ID1','stn')]

##### save the metadata of exposure and plasma outcome ####

write.table(species_filt, file = "./data/species_filt.csv", sep = ",", col.names = TRUE, row.names = TRUE)





##### metadata by stool ####


meta_stn <- all_species_filt[names(all_species_filt) == c('ID1', 'stn')]
meta_stn$id_st <- row.names(meta_stn)
meta_stn$ID1 <- as.integer(meta_stn$ID1)


##****find how to put weekly data separately ***##

w1 <- meta_pl_id %>% select(ID1, contains('_w1'))
w2 <- meta_pl_id %>% select(ID1, contains('_w2'))

colnames(w1) <- str_replace_all(colnames(w1),  c("_w1" = "_w"))
colnames(w2) <- str_replace_all(colnames(w1),  c("_w2" = "_w"))

stn12 <- meta_stn %>% filter(stn %in% c("s1", "s2"))
stn34 <- meta_stn %>% filter(stn %in% c("s3", "s4"))

id_w1 <- left_join(stn12, w1, by = "ID1")
id_w2 <- left_join(stn34, w2, by = "ID1")
id_w <- rbind(id_w1, id_w2)
id_w <- id_w %>% select(-c(ID1, stn))

no_w1w2 <- meta_pl_id %>% select(everything(), -contains('_w1'), -contains('_w2'))
meta_stn <- left_join(meta_stn, no_w1w2, by = "ID1")
meta_stn <- left_join(meta_stn, id_w, by = "id_st")
rownames(meta_stn) <- meta_stn$id_st
meta_stn <- meta_stn %>% select(-c(id_st))
# avg <- meta_id %>% select(ID1, contains('avg'), -contains('_avg'))

### what about ffq1, ffq2

write.table(meta_stn, file = "data/meta_stn.tsv", sep = "\t", col.names = TRUE, row.names = TRUE)


########## DICTIONARY #############




# read in taxonomy data

#View(ttax_rpk)


meta_species_nona<-merge(meta_stn, species_filt, by=0, rownames = TRUE)
rownames(meta_species_nona) <- meta_species_nona$Row.names

write.csv(meta_species_nona, file="./data/meta_species_nona.csv")

meta_stn_nona <- meta_species_nona %>% select(c(colnames(meta_stn)))
species_nona <- meta_species_nona %>% select(c(colnames(species_filt)))
rownames(meta_stn_nona) <- rownames(meta_species_nona)
rownames(species_nona) <- rownames(meta_species_nona)


write.csv(meta_stn_nona, file="./data/meta_stn_nona.csv")
write.csv(species_nona, file="./data/species_nona.csv")

meta_species_nona$smoke12


library(haven)
sas <- read_sas("./data/mlvs_meta_data.sas")



