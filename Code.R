# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
########################## prepare dataset ####################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
library(haven)
library(data.table)
library(dplyr)
library(tidyverse)

#setwd("/Users/jorickbater/Downloads")
getwd()
### Pre-processing steps 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
########################## (1) metadata ###############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

###################### read in channing metadata ###############################

#Exposure Data
source("noGit/Rstart.R")
z_mlvs_exposure_new <- read.table("noGit/mlvs_exposure_sunjeong.txt", 
                                  header=TRUE, sep='\t', check.names=TRUE, quote ="")
#z_mlvs_exposure_new <- merge(totom,z_mlvs_exposure_new,by="id")

z_idkey <- read.csv('noGit/idkey.csv', header=TRUE, sep = ",", stringsAsFactors = FALSE)
rownames(z_idkey) <- z_idkey$id
z_idkey <- z_idkey[order(z_idkey$id), ] 

z_mlvs_exposure <- reassignKey(z_mlvs_exposure_new)

#gsub search and replace in files -> ~~_w1avg colnames changed to ~~_w1, etc.
colnames(z_mlvs_exposure) <- gsub("w1avg","w1", colnames(z_mlvs_exposure)) 
colnames(z_mlvs_exposure) <- gsub("w2avg","w2", colnames(z_mlvs_exposure))
colnames(z_mlvs_exposure) <- gsub("plasma1","w1", colnames(z_mlvs_exposure))
colnames(z_mlvs_exposure) <- gsub("plasma2","w2", colnames(z_mlvs_exposure))

z_mlvs_agebmi <- z_mlvs_exposure [ , colnames(z_mlvs_exposure) %in% c('agemlvs', 'bmi12','act10','id')] 
#id is actually rownames


####****not working. don't have ala10v etc. as colnames****
keep_metadata_cum <- c('ala10v','epa10v','dha10v','dpa10v','calor10v', 'trans10v', 'omega610v', 'omega310v', 'omega3_noala10v')
z_metadata_cum <- z_mlvs_exposure[, keep_metadata_cum]

keep_metadata_avg <- c(
  #ffq questionnaire
  'acid_avg', 'alc_avg', 'abx_avg', 'prep_avg',
  'selfdiet_avg', 'probx_avg', 'bristolcat_avg', 'bristol_avg', 'oral_avg', 'yogurt_avg',
  #ffq energy
  'calor_avg',
  #Omega nutrients
  'ala_avg', 'epa_avg', 'dha_avg', 'dpa_avg', 'omega3_avg', 'omega3_noala_avg', 'omega6_avg','trans_avg'
)

z_metadata_avg <- z_mlvs_exposure[, keep_metadata_avg]

keep_metadata_w1 <- c(
  #questionnaire
  'acid_w1', 'alc_w1', 'abx_w1', 'prep_w1',
  'selfdiet_w1', 'probx_w1', 'bristolcat_w1', 'bristolcat5', 'bristolcat6', 'oral_w1', 'yogurt_w1',
  'bristol_w1', 'bristol5', 'bristol6',
  #ddr1
  'calor_fs_dr_w1', 
  #ddr1 nutrients
  'crp_w1','logcrp_w1', 'a_omega3_fs_dr_w1', 'a_omega3_noala_fs_dr_w1', 'loghdl_w1', 'logtg_w1', 'a_ala_fs_dr_w1', 'a_f205_fs_dr_w1', 'a_f226_fs_dr_w1', 'a_p22_5_fs_dr_w1', 'a_trn07_fo_dr_w1'
)

keep_metadata_w2 <- c(
  #questionnaire
  'acid_w2', 'alc_w2', 'abx_w2', 'prep_w2',
  'selfdiet_w2', 'probx_w2', 'bristolcat_w2', 'bristolcat7', 'bristolcat8', 'oral_w2', 'yogurt_w2',
  'bristol_w2', 'bristol7', 'bristol8',
  #ddr2 energy/fiber
  'calor_fs_dr_w2',
  #ddr2 nutrients
  'crp_w2','logcrp_w2', 'a_omega3_fs_dr_w2','a_omega3_noala_fs_dr_w2', 'loghdl_w2', 'logtg_w2', 'a_ala_fs_dr_w2', 'a_f205_fs_dr_w2', 'a_f226_fs_dr_w2', 'a_p22_5_fs_dr_w2', 'a_trn07_fo_dr_w2'
)

z_metadata_w1 <- z_mlvs_exposure[, keep_metadata_w1]
z_metadata_w2 <- z_mlvs_exposure[, keep_metadata_w2]
z_metadata_all <- cbind(z_mlvs_agebmi, z_metadata_cum, z_metadata_avg, z_metadata_w1, z_metadata_w2)
#z_metadata_all <- cbind(z_mlvs_agebmi, z_metadata_avg, z_metadata_w1, z_metadata_w2)

#### subset metadata to prep for maaslin.pcl files 

###FFQ###
demographics_ffq <- c('agemlvs', 'bmicatmlvs')
questionnaire_ffq <- c('acid_avg', 'alc_avg', 'abx_avg', 'prep_avg', 'selfdiet_avg', 'probx_avg', 'bristolcat_avg', 'oral_avg', 'yogurt_avg')
energyfiber_ffqcum <- c('ala10v','epa10v','dha10v','dpa10v','calor10v', 'trans10v', 'omega610v', 'omega310v', 'omega3_noala10v') 
energyfiber_ffq <- c('calor_avg', 'ala_avg', 'epa_avg', 'dha_avg', 'dpa_avg', 'omega3_avg', 'omega3_noala_avg', 'omega6_avg','trans_avg')

###W1###
demographics_w1 <- c('agemlvs', 'bmicatmlvs') 
questionnaire_w1 <- c('acid_w1', 'alc_w1', 'abx_w1', 'prep_w1', 'selfdiet_w1', 'probx_w1', 'bristolcat_w1', 'bristolcat5', 'bristolcat6', 'oral_w1', 'yogurt_w1')
energyfiber_w1 <- c('calor_fs_dr_w1', 'a_omega3_fs_dr_w1', 'a_omega3_noala_fs_dr_w1', 'a_ala_fs_dr_w1', 'a_f205_fs_dr_w1', 'a_f226_fs_dr_w1', 'a_p22_5_fs_dr_w1', 'a_trn07_fo_dr_w1') 

###W2###
demographics_w2 <- c('agemlvs', 'bmicatmlvs')
questionnaire_w2 <- c('acid_w2', 'alc_w2', 'abx_w2', 'prep_w2', 'selfdiet_w2', 'probx_w2', 'bristolcat_w2', 'bristolcat7', 'bristolcat8', 'oral_w2', 'yogurt_w2')
energyfiber_w2 <- c('calor_fs_dr_w2', 'a_omega3_fs_dr_w2', 'a_omega3_noala_fs_dr_w2', 'a_ala_fs_dr_w2', 'a_f205_fs_dr_w2', 'a_f226_fs_dr_w2', 'a_p22_5_fs_dr_w2', 'a_trn07_fo_dr_w2') 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############################## (2) species data ###########################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# bring in legacy file - just for checking, I actually used species data in Channing
#library(dplyr)
all_species <- read_tsv("metaphlan_taxonomic_profiles.tsv")
names(all_species)[names(all_species) == '# taxonomy'] <- 'Sample'

all_species <- all_species %>% 
  rename_at(.vars = vars(ends_with("_dna_taxonomic_profile")),
            .funs = funs(sub("[_]dna_taxonomic_profile", "", .)))

all_species <-all_species %>%
  separate(Sample, c("kingdom",       "phylum",        "class" ,        "order",         "family",        "genus" ,        "species" ), 
           sep = '\\|', remove = TRUE)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#all_species <-   read.table( "./bugs_dna_929_unFilt.tsv",
                              #sep = '\t',    header = TRUE,    check.names = FALSE,    na.strings = c("", "NA"))
#all_species<-all_species %>%
  #separate(Sample, c("kingdom",       "phylum",        "class" ,        "order",         "family",        "genus" ,        "species" ,      "strain"), 
           #sep = '\\|', remove = TRUE)

# only keep species-level features
#all_species1 <- subset(all_species,!is.na(species) & is.na(strain))
all_species1 <- subset(all_species,!is.na(species))
all_species1 <- as.data.frame(all_species1) #weird formatting change to maintain rows
rownames(all_species1)<-all_species1$species
all_species1<-all_species1[,-c(1:7)]
#all_species1<-all_species1[,-c(1:8)]

all_species <- as.data.frame(all_species1)

dim(all_species)
# [1] 468 929

all_sample <- colnames(all_species)
#View(all_sample)

# remove guy with colectomy 
all_species <- all_species[ , !colnames(all_species) %in% c(grep("15074981", colnames(all_species), value = T))]
dim(all_species)
# [1] 468 925
#colSums(all_species)

# generate list of species after abundance/prevalence filtering
# since this file is scaled to 100, i am keeping bugs with >10% prevalence and 0.01 abundance (i.e. .0001*100 when relab is scaled 0-1 and not 0-100)
dim(all_species[apply(all_species, 1, function(bug) sum(bug >= 0.01) >= 0.1*ncol(all_species)), ])
# [1] 139 925 = the 139 species after QC included in initial starr manuscripts
all_species_filt <- all_species[ apply(all_species, 1, function(bug) sum(bug >= 0.01) >= 0.1*ncol(all_species)), ]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# important! bug tables from channing need to be converted to 0-1 (not 0-11), which will wreck our downstream stats

all_species_filt   <- all_species_filt  / 100 

# let's remove all the stool suffixes from the sas names to facilitate changing the name back to species names 
#colnames(z_species.s1) <- gsub("*_s\\d","", colnames(z_species.s1))

# the species i want to keep are designated by filtering above (species_list), gonna lower case it and get rid of leading "s__"
rownames(all_species_filt) <- tolower(gsub(pattern = "s__", "", rownames(all_species_filt)))
all_species_filt <- t(all_species_filt)
all_species_filt <- as.data.frame(all_species_filt)

all_species_filt$id <- substr(rownames(all_species_filt), start = 1, stop = 6)
all_species_filt$st <- substr(rownames(all_species_filt), start = 10, stop = 13)
all_species_filt$st <- str_replace_all(all_species_filt$st,  c("SF05" = "s1","SF06" = "s2","SF07" = "s3","SF08" = "s4"))
rownames(all_species_filt) <- paste0(all_species_filt$id,"_", all_species_filt$st)
z_species.all = all_species_filt
names(z_species.all)[names(z_species.all) == 'id'] <- 'participant'
names(z_species.all)[names(z_species.all) == 'st'] <- 'stoolnum'
z_species.all$stoolnum <- str_replace_all(z_species.all$stoolnum,  c("s" = ""))
species_all <- all_species_filt[,!colnames(all_species_filt) %in% c('id','st')]


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#append metadata onto short- and long-term .pcl files
# metadata_all doesn't have a participants variable 
z_metadata_all$participant <- rownames(z_metadata_all)

# merge the species and metadata tables 
z_merged.all <- merge(x = z_species.all, y = z_metadata_all, by = 'participant', all.x = TRUE)

# rename the columns so they're unique participant_s#
rownames(z_merged.all) <- paste(z_merged.all$participant, z_merged.all$stoolnum, sep='_s') #naming convention

#checkit <- colnames(z_merged.all)
#View(checkit)


# keep only related metadata to the time course we are interested in 
# reorder so that metadata on top, bacterial data on the bottom 
z_ffq.all_species.all <- z_merged.all %>% select(participant, stoolnum, colnames(z_mlvs_agebmi), colnames(z_metadata_avg), everything(), -contains('w1'), -contains('w2'))

z_w1.all_species.all <- z_merged.all %>% select(participant, stoolnum, colnames(z_mlvs_agebmi), colnames(z_metadata_w1), contains ('w1'), everything(), -contains('avg'), -contains('w2'))
z_w2.all_species.all <- z_merged.all %>% select(participant, stoolnum, colnames(z_mlvs_agebmi), colnames(z_metadata_w2), contains ('w2'), everything(), -contains('avg'), -contains('w1'))

# subset stool
z_s1.all_species.all <- subset (z_w1.all_species.all, stoolnum == 1 )
z_s1.all_species.all <- z_s1.all_species.all %>% select(-contains('bristolcat6'), -contains('bristolcat_w1'), -contains('bristolcat7'), -contains('bristolcat8'), 
                                                        -contains('bristol6'), -contains('bristol7'), -contains('bristol8'), -contains('bristol_w1'))

z_s2.all_species.all <- subset (z_w1.all_species.all, stoolnum == 2 )
z_s2.all_species.all <- z_s2.all_species.all %>% select(-contains('bristolcat5'), -contains('bristolcat_w1'), -contains('bristolcat7'), -contains('bristolcat8'),
                                                        -contains('bristol5'), -contains('bristol7'), -contains('bristol8'), -contains('bristol_w1'))

z_s3.all_species.all <- subset (z_w2.all_species.all, stoolnum == 3 )
z_s3.all_species.all <- z_s3.all_species.all %>% select(-contains('bristolcat8'), -contains('bristolcat_w2'), -contains('bristolcat5'), -contains('bristolcat6'),
                                                        -contains('bristol5'), -contains('bristol6'), -contains('bristol8'), -contains('bristol_w2'))

z_s4.all_species.all <- subset (z_w2.all_species.all, stoolnum == 4 )
z_s4.all_species.all <- z_s4.all_species.all %>% select(-contains('bristolcat7'), -contains('bristolcat_w2'), -contains('bristolcat5'), -contains('bristolcat6'),
                                                        -contains('bristol5'), -contains('bristol6'), -contains('bristol7'), -contains('bristol_w2'))

checkit <- colnames(z_s2.all_species.all)
#View(checkit)

# finally, you've got one dataset with w1 stool and w1 metadata and another with w2 stool with w2 metadata 
# now we can combine them back into one, large flat table to do the short-term analysis but one barrier is that the metadata in each 
# has a _w1 or _w2 suffix. to fix this, rename the metadata variables to remove this so you can rbind 

colnames(z_s1.all_species.all) <- gsub("_w1","", colnames(z_s1.all_species.all))
colnames(z_s1.all_species.all) <- gsub("cat5","cat", colnames(z_s1.all_species.all))
colnames(z_s1.all_species.all) <- gsub("bristol5","bristol", colnames(z_s1.all_species.all))

colnames(z_s2.all_species.all) <- gsub("_w1","", colnames(z_s2.all_species.all))
colnames(z_s2.all_species.all) <- gsub("cat6","cat", colnames(z_s2.all_species.all))
colnames(z_s2.all_species.all) <- gsub("bristol6","bristol", colnames(z_s2.all_species.all))

colnames(z_s3.all_species.all) <- gsub("_w2","", colnames(z_s3.all_species.all))
colnames(z_s3.all_species.all) <- gsub("cat7","cat", colnames(z_s3.all_species.all))
colnames(z_s3.all_species.all) <- gsub("bristol7","bristol", colnames(z_s3.all_species.all))

colnames(z_s4.all_species.all) <- gsub("_w2","", colnames(z_s4.all_species.all))
colnames(z_s4.all_species.all) <- gsub("cat8","cat", colnames(z_s4.all_species.all))
colnames(z_s4.all_species.all) <- gsub("bristol8","bristol", colnames(z_s4.all_species.all))

# 
# z_w12.all_species.all <- rbind (z_w1.all_species.all, z_w2.all_species.all)
z_w12.all_species.all <- rbind (z_s1.all_species.all, z_s2.all_species.all, z_s3.all_species.all, z_s4.all_species.all)

checkit <- colnames(z_w12.all_species.all)
View(checkit)

# BEAUTIFUL! now you have full metadata :: species tables for both long- and short-term

# 1. z_ffq.all_species.all = long-term :: species
# 2. z_w12.all_species.all = short-term :: species 

# metadata set
# cohort FFQ cumulative average
#'ala10v','epa10v','dha10v', 'dpa10v', 'trans10v', 'omega610v', 'omega310v', 'omega3_noala10v'
# Note: I tried further including act10 and bristol score as continuous in the model, which did not influence the results too much.
ffqcum_keep <- c('participant', 'agemlvs', 'calor10v', 'abx_avg')
ffqcum_keep_nocal <- c('participant', 'agemlvs', 'abx_avg')
ala_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep, 'ala10v')] 
ala_nocal_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep_nocal, 'ala10v')] 
epa_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep, 'epa10v')] 
epa_nocal_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep_nocal, 'epa10v')] 
dha_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep, 'dha10v')] 
dha_nocal_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep_nocal, 'dha10v')] 
dpa_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep, 'dpa10v')] 
dpa_nocal_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep_nocal, 'dpa10v')] 
trans_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep, 'trans0v')] 
trans_nocal_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep_nocal, 'trans10v')] 
omega6_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep, 'omega610v')] 
omega6_nocal_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep_nocal, 'omega610v')] 
omega3_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep, 'omega310v')] 
omega3_nocal_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep_nocal, 'omega310v')] 
omega3_noala_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep, 'omega3_noala10v')] 
omega3_noala_nocal_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep_nocal, 'omega3_noala10v')] 

# MLVS FFQ
#'ala_avg', 'epa_avg', 'dha_avg', 'dpa_avg', 'omega3_avg', 'omega3_noala_avg', 'omega6_avg','trans_avg'
ffq_keep <- c('participant', 'agemlvs', 'calor_avg', 'abx_avg')
ffq_keep_nocal <- c('participant', 'agemlvs', 'abx_avg')

ala_avg_fs_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep, 'ala_avg')] 
ala_avg_fs_nocal_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep_nocal, 'ala_avg')] 

epa_avg_fs_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep, 'epa_avg')] 
epa_avg_fs_nocal_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep_nocal, 'epa_avg')] 

dha_avg_fs_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep, 'dha_avg')] 
dha_avg_fs_nocal_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep_nocal, 'dha_avg')] 

dpa_avg_fs_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep, 'dpa_avg')] 
dpa_avg_fs_nocal_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep_nocal, 'dpa_avg')] 

omega3_avg_fs_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep, 'omega3_avg')] 
omega3_avg_fs_nocal_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep_nocal, 'omega3_avg')] 

omega3_noala_avg_fs_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep, 'omega3_noala_avg')] 
omega3_noala_avg_fs_nocal_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep_nocal, 'omega3_noala_avg')] 

omega6_avg_fs_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep, 'omega6_avg')] 
omega6_avg_fs_nocal_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep_nocal, 'omega6_avg')] 

trans_avg_fs_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep, 'trans_avg')] 
trans_avg_fs_nocal_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep_nocal, 'trans_avg')] 

# MLVS DDR
#'a_omega3_fs', 'a_omega3_noala_fs', 'a_ala_fs_dr_w1', 'a_f205_fs_dr_w1', 'a_f226_fs_dr_w1', 'a_p22_5_fs_dr_w1', 'a_trn07_fo_dr_w1'
ddr_keep <- c('participant', 'agemlvs', 'calor_fs_dr', 'abx')
ddr_keep_nocal <- c('participant', 'agemlvs', 'abx')
a_omega3_fs_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_omega3_fs_dr')] 
a_omega3_fs_nocal_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep_nocal, 'a_omega3_fs_dr')] 
a_omega3_noala_fs_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_omega3_noala_fs_dr')] 
a_omega3_noala_fs_nocal_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep_nocal, 'a_omega3_noala_fs_dr')] 
a_ala_fs_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_ala_fs_dr')] 
a_ala_fs_nocal_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep_nocal, 'a_ala_fs_dr')] 
a_trans_fs_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_trn07_fo_dr')] 
a_trans_fs_nocal_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep_nocal, 'a_trn07_fo_dr')] 
a_epa_fs_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_f205_fs_dr')] 
a_epa_fs_nocal_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep_nocal, 'a_f205_fs_dr')] 
a_dha_fs_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_f226_fs_dr')] 
a_dha_fs_nocal_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep_nocal, 'a_f226_fs_dr')] 
a_dpa_fs_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_p22_5_fs_dr')] 
a_dpa_fs_nocal_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep_nocal, 'a_p22_5_fs_dr')] 
omega3_tfat_ddr <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_omega3_fs_dr', "a_fat_fs_dr")] 
trans_tfat_ddr <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_trans_fs_dr', "a_fat_fs_dr")] 
epa_tfat_ddr <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_epa_fs_dr', "a_fat_fs_dr")] 
dha_tfat_ddr <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_dha_fs_dr', "a_fat_fs_dr")] 
dpa_tfat_ddr <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_dpa_fs_dr', "a_fat_fs_dr")] 

# MLVS biomarker
biomarker_keep <- c('participant', 'agemlvs', 'abx', 'bmi12')
biomarker_logcrp  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logcrp')] 
biomarker_logcrp_complete <- biomarker_logcrp[!is.na(biomarker_logcrp$logcrp), ]
biomarker_loghdl  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'loghdl')] 
biomarker_loghdl_complete <- biomarker_loghdl[!is.na(biomarker_loghdl$loghdl), ]
biomarker_logtg  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtg')] 
biomarker_logtg_complete <- biomarker_logtg[!is.na(biomarker_logtg$logtg), ]
a_omega3_fs_ddr_logcrp <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logcrp','a_omega3_fs_dr')] 
a_omega3_fs_ddr_logtchdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtchdl','a_omega3_fs_dr')] 
a_omega3_fs_ddr_loghdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'loghdl','a_omega3_fs_dr')] 
a_omega3_fs_ddr_logtg <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtg','a_omega3_fs_dr')] 

a_omega3_noala_fs_ddr_logcrp <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logcrp','a_omega3_noala_fs_dr')] 
a_omega3_noala_fs_ddr_logtchdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtchdl','a_omega3_noala_fs_dr')] 
a_omega3_noala_fs_ddr_loghdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'loghdl','a_omega3_noala_fs_dr')] 
a_omega3_noala_fs_ddr_logtg <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtg','a_omega3_noala_fs_dr')] 

a_ala_fs_ddr_logcrp <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logcrp','a_ala_fs_dr')] 
a_ala_fs_ddr_logtchdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtchdl','a_ala_fs_dr')] 
a_ala_fs_ddr_loghdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'loghdl','a_ala_fs_dr')] 
a_ala_fs_ddr_logtg <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtg','a_ala_fs_dr')] 

a_trans_fs_ddr_logcrp <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logcrp','a_trans_fs_dr')] 
a_trans_fs_ddr_logtchdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtchdl','a_trans_fs_dr')] 
a_trans_fs_ddr_loghdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'loghdl','a_trans_fs_dr')] 
a_trans_fs_ddr_logtg <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtg','a_trans_fs_dr')] 

a_epa_fs_ddr_logcrp <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logcrp','a_epa_fs_dr')] 
a_epa_fs_ddr_logtchdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtchdl','a_epa_fs_dr')] 
a_epa_fs_ddr_loghdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'loghdl','a_epa_fs_dr')] 
a_epa_fs_ddr_logtg <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtg','a_epa_fs_dr')] 

a_dha_fs_ddr_logcrp <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logcrp','a_dha_fs_dr')] 
a_dha_fs_ddr_logtchdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtchdl','a_dha_fs_dr')] 
a_dha_fs_ddr_loghdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'loghdl','a_dha_fs_dr')] 
a_dha_fs_ddr_logtg <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtg','a_dha_fs_dr')] 

a_dpa_fs_ddr_logcrp <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logcrp','a_dpa_fs_dr')] 
a_dpa_fs_ddr_logtchdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtchdl','a_dpa_fs_dr')] 
a_dpa_fs_ddr_loghdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'loghdl','a_dpa_fs_dr')] 
a_dpa_fs_ddr_logtg <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtg','a_dpa_fs_dr')] 

#drop observations with missing short-term fiber
#I don't know why I need to use the complete dataset to run MaAsLin2 in my mac...
a_omega3_fs_ddr_complete <- a_omega3_fs_ddr[a_omega3_fs_ddr$a_omega3_fs_dr>0,]
a_omega3_fs_ddr_complete <- na.omit(a_omega3_fs_ddr_complete)

# Tobyn
# check distribution of fiber 
graphics.off()
histogram(a_omega3_fs_ddr$a_omega3_fs_dr)
quantile(a_omega3_fs_ddr$a_omega3_fs_dr,  probs = seq (0,1,0.1), na.rm = TRUE,
         names = TRUE)
#       0%       10%       20%       30%       40%       50%       60%       70% 
#9.392454 16.591799 19.089333 20.941394 22.565379 24.031779 25.967835 27.866991 
#80%       90%      100% 
#30.168380 35.819875 55.518952 

IQR(a_omega3_fs_ddr$a_omega3_fs_dr, na.rm = TRUE) 
#1.085002
quantile(a_omega3_fs_ddr$a_omega3_fs_dr,  0.75, na.rm = TRUE,
         names = TRUE)
#     75% 
# 2.75162
# outlier: Q3+1.5IQR=2.75162+1.5*1.085002=4.379123
# Q3+3IQR=2.75162+3*1.085002=6.006626

#drop observations with short-term fiber >Q3+1.5IQR
a_omega3_fs_ddr_1.5IQR <- a_omega3_fs_ddr[a_omega3_fs_ddr$a_omega3_fs_dr<=4.379123,]
a_omega3_fs_ddr_3IQR <- a_omega3_fs_ddr[a_omega3_fs_ddr$a_omega3_fs_fs_dr<=6.006626,]

mlvs_metadata_925 <- z_w12.all_species.all[,!colnames(z_w12.all_species.all) %in% colnames(species_all)]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#############  pcl writing  ############
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

for(i in 1:ncol(species_all)) 
{ 
  species_all[ , i] <- as.numeric(as.character(species_all[, i])) 
}

pcl_file <- c('species_all', 
              'ala_ffqcum','epa_ffqcum','dha_ffqcum', 'dpa_ffqcum', 'trans_ffqcum', 'omega6_ffqcum', 'omega3_ffqcum', 'omega3_noala_ffqcum',
              'ala_avg_fs_ffq', 'epa_avg_fs_ffq', 'dha_avg_fs_ffq', 'dpa_avg_fs_ffq', 'trans_avg_fs_ffq', 'omega6_avg_fs_ffq', 'omega3_avg_fs_ffq', 'omega3_noala_avg_fs_ffq',
              'a_omega3_fs_ddr', 'a_omega3_fs_ddr_complete', 'a_omega3_fs_ddr_1.5IQR', 'a_omega3_fs_ddr_3IQR', 'a_omega3_noala_fs_ddr', 'a_ala_fs_ddr', 'a_trans_fs_ddr', 'a_epa_fs_ddr', 'a_dha_fs_ddr', 'a_dpa_fs_ddr',
              'omega3_tfat_ddr', 'trans_tfat_ddr', 'epa_tfat_ddr', 'dha_tfat_ddr', 'dpa_tfat_ddr',
              'biomarker_logcrp', 'biomarker_logcrp_complete', 'a_omega3_fs_ddr_logcrp', 'biomarker_loghdl', 'biomarker_loghdl_complete', 'a_omega3_fs_ddr_loghdl', 'biomarker_logtg', 'biomarker_logtg_complete', 'a_omega3_fs_ddr_logtg',
              'a_omega3_noala_fs_ddr_logcrp', 'a_omega3_noala_fs_ddr_loghdl', 'a_omega3_noala_fs_ddr_logtg',
              'a_trans_fs_ddr_logcrp', 'a_trans_fs_ddr_loghdl', 'a_trans_fs_ddr_logtg',
              'a_ala_fs_ddr_logcrp', 'a_ala_fs_ddr_loghdl', 'a_ala_fs_ddr_logtg',
              'a_epa_fs_ddr_logcrp', 'a_epa_fs_ddr_loghdl', 'a_epa_fs_ddr_logtg',
              'a_dha_fs_ddr_logcrp', 'a_dha_fs_ddr_loghdl', 'a_dha_fs_ddr_logtg',
              'a_dpa_fs_ddr_logcrp', 'a_dpa_fs_ddr_loghdl', 'a_dpa_fs_ddr_logtg',
              'z_metadata_all','mlvs_metadata_925')   

# make pcl files 
for (a in 1:length(pcl_file))
{
  tempdf <- NULL
  # transpose
  tempdf <- as.data.frame ( t (get( pcl_file[a] ) ) )
  # steps to add the word sample to cell [1,1] so maaslin won't throw an error 
  tempcol <- colnames(tempdf)
  tempdf$sample = row.names(tempdf)
  addsamp = c('sample', tempcol)
  foo = tempdf [ , addsamp]
  
  
  
  # write table
  write.table(foo, paste('./maaslin2/pcl/', pcl_file[a], '.pcl',sep=''),  sep = '\t', quote = F, eol = '\n',  row.names=F)
}   

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
source("/Users/jorickbater/Desktop/Mingyang Stuff/Data/Rstart.R")

# read in metadata (925)
mlvs_metadata_925 <- read.table('./maaslin2/pcl/mlvs_metadata_925.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(mlvs_metadata_925) <- mlvs_metadata_925 $sample
mlvs_metadata_925  <- mlvs_metadata_925 [,-1]
mlvs_metadata_925  <- as.data.frame(t(mlvs_metadata_925 ))

# read in species
species_all <- read.table('./maaslin2/pcl/species_all.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(species_all) <- species_all$sample
species_all <- species_all[,-1]
species_all <- as.data.frame(t(species_all))

# merge metadata and species
species_metadata <- merge(x = mlvs_metadata_925, y = species_all, by = "row.names", all = TRUE)
rownames(species_metadata) <- species_metadata[ ,1]
species_metadata[ ,1] <- NULL

ddr_plot <- species_metadata
checkit <- colnames(ddr_plot)
#View(checkit)

for(i in 1:ncol(ddr_plot)) 
{ 
  ddr_plot[ , i] <- as.numeric(as.character(ddr_plot[, i])) 
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
####################### Fig 3A. MaAsLin2 ################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

bug_file<-c('species_all')
metadata_file <- c( 'ala_ffqcum','epa_ffqcum','dha_ffqcum', 'dpa_ffqcum', 'trans_ffqcum', 'omega6_ffqcum', 'omega3_ffqcum', 'omega3_noala_ffqcum',
                    'ala_avg_fs_ffq', 'epa_avg_fs_ffq', 'dha_avg_fs_ffq', 'dpa_avg_fs_ffq', 'trans_avg_fs_ffq', 'omega6_avg_fs_ffq', 'omega3_avg_fs_ffq', 'omega3_noala_avg_fs_ffq', 'a_omega3_fs_ddr', 'a_omega3_noala_fs_ddr', 'a_trans_fs_ddr', 'a_ala_fs_ddr', 'a_epa_fs_ddr', 'a_dha_fs_ddr', 'a_dpa_fs_ddr',
                   'biomarker_logcrp', 'biomarker_logcrp_complete', 'a_omega3_fs_ddr_logcrp', 'biomarker_loghdl', 'biomarker_loghdl_complete', 'a_omega3_fs_ddr_loghdl', 'biomarker_logtg', 'biomarker_logtg_complete', 'a_omega3_fs_ddr_logtg', 'a_trans_fs_ddr_logcrp', 'a_trans_fs_ddr_loghdl', 'a_trans_fs_ddr_logtg',
                   'a_ala_fs_ddr_logcrp', 'a_ala_fs_ddr_loghdl', 'a_ala_fs_ddr_logtg',
                   'a_epa_fs_ddr_logcrp', 'a_epa_fs_ddr_loghdl', 'a_epa_fs_ddr_logtg',
                   'a_dha_fs_ddr_logcrp', 'a_dha_fs_ddr_loghdl', 'a_dha_fs_ddr_logtg',
                   'a_dpa_fs_ddr_logcrp', 'a_dpa_fs_ddr_loghdl', 'a_dpa_fs_ddr_logtg')


# arcsin sqrt transformation  
for (a in 1:length(metadata_file))
{
  # (1) with no abundance/prevalence filter since i already did this myself
  
  for (b in 1:length(bug_file))
  {
    Maaslin2(input_data     = paste('./maaslin2/pcl/', bug_file[b], '.pcl', sep=''), 
             input_metadata   = paste('./maaslin2/pcl/', metadata_file[a], '.pcl', sep=''),
             output           = paste('./maaslin2/omega output', metadata_file[a], '_', bug_file[b], '/', sep=''),
             normalization    = 'NONE', 
             standardize      = 'FALSE',
             transform        = 'AST', 
             analysis_method  = 'LM', 
             random_effects   = 'participant',
             min_abundance    = 0, 
             min_prevalence   = 0,
             cores            = 1)
    
  }
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################## Figure 3B/C. scatter plots ###################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function of positive associations
scatter_plot_positive <- function(x,y,xlab,ylab) {
  plot <- ggplot(
    data=ddr_plot, aes_string(x=x, y=y)) +
    geom_point( aes(), fill = '#D9EF8B', color = 'black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
    stat_smooth(method = "glm", color ='black')+ 
    guides(alpha='none')+labs("")+
    xlab(xlab) +  ylab(ylab)  +
    nature_theme +
    guides(legend.position=NULL)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title=element_text(size=20),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.ticks =element_blank(),
          plot.title = element_blank())
  filepath <- paste('./result/figures/scatter_species/', y, '_', x, '.png', sep ="" )
  ggsave(filename=filepath, plot=plot, width = 5, height = 5, dpi = 600) 
}# end of function

scatter_plot_positive('a_aofib_fs_dr','eubacterium_eligens',"Short-term dietary fiber\n(g/d)","Eubacterium eligens\n(relative abundance")
scatter_plot_positive('a_aofib_fs_dr','faecalibacterium_prausnitzii',"Short-term dietary fiber\n(g/d)","Faecalibacterium prausnitzii\n(relative abundance")
scatter_plot_positive('logcrp','bacteroides_uniformis',"Log CRP","Bacteroides uniformis\n(relative abundance)")


# Function of negative associations
scatter_plot_negative <- function(x,y,xlab,ylab) {
  plot <- ggplot(
    data=ddr_plot, aes_string(x=x, y=y)) +
    geom_point( aes(), fill = '#C92D39', color = 'black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
    stat_smooth(method = "glm", color ='black')+ 
    guides(alpha='none')+labs("")+
    xlab(xlab) +  ylab(ylab)  +
    nature_theme +
    guides(legend.position=NULL)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title=element_text(size=20),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.ticks =element_blank(),
          plot.title = element_blank())
  filepath <- paste('./result/figures/scatter_species/', y, '_', x, '.png', sep ="" )
  ggsave(filename=filepath, plot=plot, width = 5, height = 5, dpi = 600) 
}# end of function

scatter_plot_negative('a_aofib_fs_dr','lachnospiraceae_bacterium_1_4_56faa',"Short-term dietary fiber\n(g/d)","Lachnospiraceae bacterium 1 4 56FAA\n(relative abundance")
scatter_plot_negative('a_aofib_fs_dr','ruminococcus_torques',"Short-term dietary fiber\n(g/d)","Ruminococcus torques\n(relative abundance")
scatter_plot_negative('logcrp','eubacterium_eligens',"Short-term dietary fiber\n(g/d)","Eubacterium eligens\n(relative abundance")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
######### Figure 2. interaction between DDR fiber and P. copri on logCRP#####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################# read in data ###################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

source("/Users/jorickbater/Desktop/Mingyang Stuff/Data/Rstart.R")

# read in metadata (925)
mlvs_metadata_925 <- read.table('./maaslin2/pcl/mlvs_metadata_925.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(mlvs_metadata_925) <- mlvs_metadata_925 $sample
mlvs_metadata_925  <- mlvs_metadata_925 [,-1]
mlvs_metadata_925  <- as.data.frame(t(mlvs_metadata_925 ))

# read in species
species_all <- read.table('./maaslin2/pcl/species_all.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(species_all) <- species_all$sample
species_all <- species_all[,-1]
species_all <- as.data.frame(t(species_all))

# read in MDS1 and MDS2
bugs_pcoa <- species_all
bugs_pcoa <- capscale(t(bugs_pcoa) ~ 1, distance = 'bray') #Question? old pcoa used asin sqrt transformed bug data?
variance.bugs <- head(eigenvals(bugs_pcoa)/sum(eigenvals(bugs_pcoa)))
x_variance.bugs <- as.integer(variance.bugs[1]*100)
y_variance.bugs <- as.integer(variance.bugs[2]*100)

# sample scores
ord.scores <- as.data.frame(scores(bugs_pcoa, display = "sites"))

# bug scores
ord.bug.scores <- as.data.frame(scores(bugs_pcoa, display = "species") )
#View(abs(ord.bug.scores))

# merge the data (change names for trans)
association <- mlvs_metadata_925[ , colnames(mlvs_metadata_925) %in% c('participant', 'agemlvs','abx', 'bmi12', 'logcrp','loghdl','logtg','ala10v','epa10v','dha10v', 'dpa10v', 'calor10v', 'a_omega3_noala_fs_dr', 'calor_fs_dr','act10', 'bristolcat')] 

#merge association data with MDS bug scores
association_mds <- merge(x = association, y = ord.bug.scores, by = "row.names", all = TRUE)

# set first column as rownames and remove it
rownames(association_mds) <- association_mds[ ,1]
association_mds[ ,1] <- NULL

# merge with the bug data
mns <- colMeans(species_all, na.rm=TRUE)
order(-mns)
species_all_sort <- species_all[,order(-mns)]
colMeans(species_all_sort)

association_bugs <- merge(x = association_mds, y = species_all_sort, by = 'row.names', all = TRUE)
rownames(association_bugs) <- association_bugs[ ,1]
association_bugs[ ,1] <- NULL

head(association_bugs)
#sort by MDS1
association_bugs_sort <- association_bugs[order(association_bugs$MDS1),]

association_bugs$prevotella_copri2c <- with(association_bugs, ifelse(prevotella_copri>0, 1, 0))
table(association_bugs$prevotella_copri2c)
#0   1 
#693 232

#association_bugs$lactobacillus_rogosae <- with(association_bugs, ifelse(lactobacillus_rogosae>0, 1, 0))
#association_bugs$bifidobacterium_adolescentis <- with(association_bugs, ifelse(bifidobacterium_adolescentis>0, 1, 0))
#association_bugs$bifidobacterium_longum<- with(association_bugs, ifelse(bifidobacterium_longum>0, 1, 0))
#association_bugs$bifidobacterium_pseudocatenulatum<-- with(association_bugs, ifelse(bifidobacterium_pseudocatenulatum>0, 1, 0))
#association_bugs$bifidobacterium_bifidum<- with(association_bugs, ifelse(bifidobacterium_bifidum>0, 1, 0))

library(lmerTest)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### interaction model with MDS1 ####################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#crude interaction model: n=775, p-int=0.663 (1), 0.031 (2)
lm.crp <- lmer(formula = logcrp ~ a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + (1 | participant), data=association_bugs)
anova(lm.crp)
summary(lm.crp)

lm.crp <- lmer(formula = logcrp ~ a_omega3_fs_dr + MDS2 + a_omega3_fs_dr*MDS2 + (1 | participant), data=association_bugs)
anova(lm.crp)
summary(lm.crp)

#further adjusting for age: n=775, p-int=0.689 (1), 0.037 (2)
lm.crp <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + (1 | participant), data=association_bugs)
anova(lm.crp)
summary(lm.crp)

lm.crp <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + MDS2 + a_omega3_fs_dr*MDS2 + (1 | participant), data=association_bugs)
anova(lm.crp)
summary(lm.crp)

#main model, try further adjusting for some covariates: age, antibiotics, calorie; n=773, p-int=0.749 (1), 0.044 (2)
lm.crp <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + calor_fs_dr + abx + (1 | participant), data=association_bugs)
anova(lm.crp)
summary(lm.crp)

lm.crp <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + MDS2 + a_omega3_fs_dr*MDS2 + calor_fs_dr + abx + (1 | participant), data=association_bugs)
anova(lm.crp)
summary(lm.crp)

#additionally adjusting for bristolcat or physical activity did not change the results; n=752, p-int=0.9113 (1), 0.085 (2)
association_bugs$bristolcat <- as.factor(association_bugs$bristolcat)
lm.crp <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + calor_fs_dr + abx + bristolcat + act10 + (1 | participant), data=association_bugs)
anova(lm.crp)
summary(lm.crp)

lm.crp <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + MDS2 + a_omega3_fs_dr*MDS2 + calor_fs_dr + abx + bristolcat + act10 + (1 | participant), data=association_bugs)
anova(lm.crp)
summary(lm.crp)

############HDL
#crude interaction model: n=775, p-int=0.663
lm.hdl <- lmer(formula = loghdl ~ a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + (1 | participant), data=association_bugs)
anova(lm.hdl)
summary(lm.hdl)

#further adjusting for age: n=775, p-int=0.689
lm.hdl <- lmer(formula = loghdl ~ agemlvs + a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + (1 | participant), data=association_bugs)
anova(lm.hdl)
summary(lm.hdl)

#main model, try further adjusting for some covariates: age, antibiotics, calorie; n=773, p-int=0.749
lm.hdl <- lmer(formula = loghdl ~ agemlvs + a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + calor_fs_dr + abx + (1 | participant), data=association_bugs)
anova(lm.hdl)
summary(lm.hdl)

#additionally adjusting for bristolcat or physical activity did not change the results; n=752, p-int=0.9113
association_bugs$bristolcat <- as.factor(association_bugs$bristolcat)
lm.hdl <- lmer(formula = loghdl ~ agemlvs + a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + calor_fs_dr + abx + bristolcat + act10 + (1 | participant), data=association_bugs)
anova(lm.hdl)
summary(lm.hdl)

############TG############
#crude interaction model: n=775, p-int=0.663
lm.tg <- lmer(formula = logtg ~ a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + (1 | participant), data=association_bugs)
anova(lm.tg)
summary(lm.tg)

#further adjusting for age: n=775, p-int=0.689
lm.tg <- lmer(formula = logtg ~ agemlvs + a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + (1 | participant), data=association_bugs)
anova(lm.tg)
summary(lm.tg)

#main model, try further adjusting for some covariates: age, antibiotics, calorie; n=773, p-int=0.749
lm.tg <- lmer(formula = logtg ~ agemlvs + a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + calor_fs_dr + abx + (1 | participant), data=association_bugs)
anova(lm.tg)
summary(lm.tg)

#additionally adjusting for bristolcat or physical activity did not change the results; n=752, p-int=0.9113
association_bugs$bristolcat <- as.factor(association_bugs$bristolcat)
lm.tg <- lmer(formula = logtg ~ agemlvs + a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + calor_fs_dr + abx + bristolcat + act10 + (1 | participant), data=association_bugs)
anova(lm.tg)
summary(lm.tg)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
########## interaction model with P copri ############
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#model adjusting for age only, p-int=0.49 ########
lm.crp.prevotella_copri <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + prevotella_copri2c + a_omega3_fs_dr*prevotella_copri2c + (1 | participant), data=association_bugs)
anova(lm.crp.prevotella_copri)
summary(lm.crp.prevotella_copri)

#lm.crp.lactobacillus_rogosae <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + lactobacillus_rogosae + a_omega3_fs_dr*lactobacillus_rogosae + (1 | participant), data=association_bugs)
#anova(lm.crp.lactobacillus_rogosae)
#summary(lm.crp.lactobacillus_rogosae)

#lm.crp.bifidobacterium_adolescentis <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + bifidobacterium_adolescentis + a_omega3_fs_dr*bifidobacterium_adolescentis + (1 | participant), data=association_bugs)
#anova(lm.crp.bifidobacterium_adolescentis)
#summary(lm.crp.bifidobacterium_adolescentis)

#lm.crp.bifidobacterium_longum <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + bifidobacterium_longum + a_omega3_fs_dr*bifidobacterium_longum + (1 | participant), data=association_bugs)
#anova(lm.crp.bifidobacterium_longum)
#summary(lm.crp.bifidobacterium_longum)

#lm.crp.bifidobacterium_pseudocatenulatum <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + bifidobacterium_pseudocatenulatum + a_omega3_fs_dr*bifidobacterium_pseudocatenulatum + (1 | participant), data=association_bugs)
#anova(lm.crp.bifidobacterium_pseudocatenulatum)
#summary(lm.crp.bifidobacterium_pseudocatenulatum)

#lm.crp.bifidobacterium_bifidum <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + bifidobacterium_bifidum + a_omega3_fs_dr*bifidobacterium_bifidum + (1 | participant), data=association_bugs)
#anova(lm.crp.bifidobacterium_bifidum)
#summary(lm.crp.bifidobacterium_bifidum)

#model adjusting for other covariates, p-int=0.46 ########
lm.crp.prevotella_copri <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + prevotella_copri2c + a_omega3_fs_dr*prevotella_copri2c + calor_fs_dr + abx + (1 | participant), data=association_bugs)
anova(lm.crp.prevotella_copri)
summary(lm.crp.prevotella_copri)

#further adjusting for bristol score and physical activity, p-int=0.30 #########
lm.crp.prevotella_copri <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + prevotella_copri2c + a_omega3_fs_dr*prevotella_copri2c + calor_fs_dr + abx + bristolcat + act10 + (1 | participant), data=association_bugs)
anova(lm.crp.prevotella_copri)
summary(lm.crp.prevotella_copri)

############HDL
#model adjusting for age only, p-int=0.33 ########
lm.hdl.prevotella_copri <- lmer(formula = loghdl ~ agemlvs + a_omega3_fs_dr + prevotella_copri2c + a_omega3_fs_dr*prevotella_copri2c + (1 | participant), data=association_bugs)
anova(lm.hdl.prevotella_copri)
summary(lm.hdl.prevotella_copri)

#model adjusting for other covariates, p-int=0.26 ########
lm.hdl.prevotella_copri <- lmer(formula = loghdl ~ agemlvs + a_omega3_fs_dr + prevotella_copri2c + a_omega3_fs_dr*prevotella_copri2c + calor_fs_dr + abx + (1 | participant), data=association_bugs)
anova(lm.hdl.prevotella_copri)
summary(lm.hdl.prevotella_copri)

#further adjusting for bristol score and physical activity, p-int=0.37 #########
lm.hdl.prevotella_copri <- lmer(formula = loghdl ~ agemlvs + a_omega3_fs_dr + prevotella_copri2c + a_omega3_fs_dr*prevotella_copri2c + calor_fs_dr + abx + bristolcat + act10 + (1 | participant), data=association_bugs)
anova(lm.hdl.prevotella_copri)
summary(lm.hdl.prevotella_copri)

############TG
#model adjusting for age only, p-int=0.99 ########
lm.tg.prevotella_copri <- lmer(formula = logtg ~ agemlvs + a_omega3_fs_dr + prevotella_copri2c + a_omega3_fs_dr*prevotella_copri2c + (1 | participant), data=association_bugs)
anova(lm.tg.prevotella_copri)
summary(lm.tg.prevotella_copri)

#model adjusting for other covariates, p-int=0.97 ########
lm.tg.prevotella_copri <- lmer(formula = logtg ~ agemlvs + a_omega3_fs_dr + prevotella_copri2c + a_omega3_fs_dr*prevotella_copri2c + calor_fs_dr + abx + (1 | participant), data=association_bugs)
anova(lm.tg.prevotella_copri)
summary(lm.tg.prevotella_copri)

#further adjusting for bristol score and physical activity, p-int=0.93 #########
lm.tg.prevotella_copri <- lmer(formula = logtg ~ agemlvs + a_omega3_fs_dr + prevotella_copri2c + a_omega3_fs_dr*prevotella_copri2c + calor_fs_dr + abx + bristolcat + act10 + (1 | participant), data=association_bugs)
anova(lm.tg.prevotella_copri)
summary(lm.tg.prevotella_copri)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############### Scatter plots by P copri ####################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#scatter plot with smooth line 
association_bugs$prevotella_copri2c <- as.factor(association_bugs$prevotella_copri2c)
head(association_bugs)

ggp_by_prevotella_copri2c <- ggplot(association_bugs, aes(x=a_omega3_fs_dr, y=logcrp, color=prevotella_copri2c))+ 
  geom_point(color="black", alpha = .5, shape = 21, size = 3, stroke = 1, na.rm=TRUE) +
  geom_point(aes(color=prevotella_copri2c), size=3)+
  geom_smooth(method=glm, se=FALSE, fullrange=TRUE)+
  labs(x="Dietary omega 3 intake (g/day)", y = "Log-transformed CRP")+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_text(size=20),
        legend.text = element_text(size=20), 
        axis.title=element_text(size=20,face='bold'),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.ticks =element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_blank())+
  scale_color_discrete(name ="Presence of Prevotella copri",labels=c("No (n=588)", "Yes (n=185)"))+
  annotate("text", x = 10, y = 2, label = "P-interaction=0.46", size=8)

plot(ggp_by_prevotella_copri2c)

############HDL
ggp_by_prevotella_copri2c <- ggplot(association_bugs, aes(x=a_omega3_fs_dr, y=loghdl, color=prevotella_copri2c))+ 
  geom_point(color="black", alpha = .5, shape = 21, size = 3, stroke = 1, na.rm=TRUE) +
  geom_point(aes(color=prevotella_copri2c), size=3)+
  geom_smooth(method=glm, se=FALSE, fullrange=TRUE)+
  labs(x="Dietary omega 3 intake (g/day)", y = "Log-transformed hdl")+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_text(size=20),
        legend.text = element_text(size=20), 
        axis.title=element_text(size=20,face='bold'),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.ticks =element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_blank())+
  scale_color_discrete(name ="Presence of Prevotella copri",labels=c("No (n=588)", "Yes (n=185)"))+
  annotate("text", x = 10, y = 2, label = "P-interaction=0.46", size=8)

plot(ggp_by_prevotella_copri2c)

############TG
ggp_by_prevotella_copri2c <- ggplot(association_bugs, aes(x=a_omega3_fs_dr, y=logtg, color=prevotella_copri2c))+ 
  geom_point(color="black", alpha = .5, shape = 21, size = 3, stroke = 1, na.rm=TRUE) +
  geom_point(aes(color=prevotella_copri2c), size=3)+
  geom_smooth(method=glm, se=FALSE, fullrange=TRUE)+
  labs(x="Dietary omega 3 intake (g/day)", y = "Log-transformed tg")+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_text(size=20),
        legend.text = element_text(size=20), 
        axis.title=element_text(size=20,face='bold'),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.ticks =element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_blank())+
  scale_color_discrete(name ="Presence of Prevotella copri",labels=c("No (n=588)", "Yes (n=185)"))+
  annotate("text", x = 10, y = 2, label = "P-interaction=0.46", size=8)

plot(ggp_by_prevotella_copri2c)

ggsave(filename='./ddr.crp.by.prevotella_copri2c.png', plot=ggp_by_prevotella_copri2c, width = 10, height = 5, dpi = 600) 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
######## Figure 1. Distribution of data and PcoA on species-level taxonomy#######
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
######### read in data########
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# read in metadata (307)
z_metadata_all <- read.table('./maaslin2/pcl/z_metadata_all.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(z_metadata_all) <- z_metadata_all$sample
z_metadata_all <- z_metadata_all[,-1]
meta_dis <- as.data.frame(t(z_metadata_all))

# read in species
species_all <- read.table('./maaslin2/pcl/species_all.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(species_all) <- species_all$sample
species_all <- species_all[,-1]
species_all <- as.data.frame(t(species_all))

# read in metadata (925)
mlvs_metadata_925 <- read.table('./maaslin2/pcl/mlvs_metadata_925.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(mlvs_metadata_925) <- mlvs_metadata_925 $sample
mlvs_metadata_925  <- mlvs_metadata_925 [,-1]
mlvs_metadata_925  <- as.data.frame(t(mlvs_metadata_925 ))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
####### Fig 1B. Distribution of DDR omega3 and CRP#######
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
graphics.off()

####distribution of ddr omega3, count########## Why week1?
aofib_ddr_distribution<-ggplot(meta_dis, aes(x=trans_avg))+ #a_omega3_fs_dr_w1
  geom_histogram(color="black", fill="grey", binwidth=1)+
  labs(x="Dietary trans fat intake (g/day)", y = "Count")+
  xlim(0,5)+ #16
  theme(axis.title=element_text(size=50,face = 'bold'),
        axis.text.x = element_text(size=50),
        axis.text.y = element_text(size=50),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid')
  )
plot(aofib_ddr_distribution)
#ggsave(filename='./result/figures/distribution/aofib_ddr_histogram_count.png', plot=aofib_ddr_distribution, width = 12, height = 5, dpi = 600)

####distribution of logcrp, count##########
summary(meta_dis$logcrp_w1)
logcrp_distribution<-ggplot(meta_dis, aes(x=logcrp_w1))+
  geom_histogram(color="black", fill="grey", binwidth=0.5)+
  labs(x="log C-reactive protein", y = "Count")+
  xlim(-2,3)+
  theme(axis.title=element_text(size=50, face = 'bold'),
        axis.text.x = element_text(size=50),
        axis.text.y = element_text(size=50),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid')
  )
plot(logcrp_distribution)
#ggsave(filename='./result/figures/distribution/logcrp_histogram_count.png', plot=logcrp_distribution, width = 12, height = 5, dpi = 600)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
######### Fig 1C. Correlation between omega3 and age/BMI##########
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

####scatter plots of ddr omega3 and age####
ddromega3_age <- ggplot(
  data=meta_dis, aes_string(x='agemlvs', y='a_omega3_fs_dr_w1')) +
  geom_point( aes(), fill = '#C92D39', color='black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
  stat_smooth(method = "glm", color ='black')+ 
  guides(alpha='none')+labs("")+
  xlab('Age (years)') +  ylab("Dietary omega3 intake (g/day)")  +
  nature_theme +
  guides(legend.position=NULL)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(size=30, face = 'bold'),
        axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=30),
        axis.ticks =element_blank(),
        plot.title = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid')
  )+
  annotate("text", x = 73, y=50, label = "spearman r=-0.02", size=15)
plot(ddromega3_age)
#ggsave(filename='./result/figures/distribution/ddrfiber_age.png', plot=ddrfiber_age, width = 10, height = 10, dpi = 600) 

cor.test(meta_dis$a_omega3_fs_dr_w1, meta_dis$agemlvs,
         alternative = c("two.sided"),
         method = c("spearman"), exact=F)



####scatter plots of ddr fiber and BMI####
ddrfiber_bmi <- ggplot(
  data=meta_dis, aes_string(x='bmi12', y='a_omega3_fs_dr_w1')) +
  geom_point( aes(), fill = '#C92D39', color='black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
  stat_smooth(method = "glm", color ='black')+ 
  guides(alpha='none')+labs("")+
  #xlab(expression(paste("Body mass index (kg/m"^"2"*")"))) +  #ylab("Dietary fiber intake (g/day)")  +
  nature_theme +
  guides(legend.position=NULL)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        #axis.title=element_text(size=30),
        #axis.title.x=element_text(size=30, face = 'bold'),
        #axis.title.y = element_blank(),
        axis.title = element_blank(),
        #axis.text.x = element_text(size=30),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=30),
        axis.ticks =element_blank(),
        plot.title = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid')
  )+
  annotate("text", x = 30, y=50, label = "spearman r=-0.24", size=15)
plot(ddrfiber_bmi)
#ggsave(filename='./result/figures/distribution/ddrfiber_bmi.png', plot=ddrfiber_bmi, width = 10, height = 10, dpi = 600) 

cor.test(meta_dis$a_omega3_fs_dr_w1, meta_dis$bmi12,
         alternative = c("two.sided"),
         method = c("spearman"), exact=F)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
######## Fig 1D. stacked barplots of fiber sources###########
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
barplot_cum <- meta_dis[,c('ala_avg','epa_avg','dha_avg','dpa_avg')]
barplot_cum <- barplot_cum %>% 
  mutate(other = ala_avg - epa_avg - dha_avg - dpa_avg)
rownames(barplot_cum) <- rownames(z_metadata_cum)
colnames(barplot_cum) <- c('total', 'ala','epa','dha','dpa')

# sort by aofib10v
barplot_cum_sort <- barplot_cum[order(-barplot_cum$total),]

# keep only fiber subsets
barplot_cum_sort_stratified <- barplot_cum_sort[,2:5]

# transform to relative abundances
#barplot_cum_sort_stratified_a <- sweep(barplot_cum_sort_stratified, 1, rowSums(barplot_cum_sort_stratified), `/`)

# transform to samples in columns and metadata in rows
t_barplot_cum_sort_stratified <- as.data.frame(t(barplot_cum_sort_stratified))


t_barplot_cum_sort_stratified$sources <- row.names(t_barplot_cum_sort_stratified)
dim(t_barplot_cum_sort_stratified) # 4 308
# melt
t_barplot_cum_sort_stratified_melt <- melt(t_barplot_cum_sort_stratified, id.vars = 'sources')
dim(t_barplot_cum_sort_stratified_melt) # 1228 3


top30_s_colors <- scale_fill_manual(values = c("#3288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D43", "#E7298A", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142", "#D53E4F", "#FDAE61", "#FEE08B", "#333333", "#66C2A5", "#3288BD", "#1B9E77", "#D95F02", "#7570B3", "#E6AB02", "#A6761D", "#CCCCCC"))

bar_cum<-ggplot(data=t_barplot_cum_sort_stratified_melt, aes(x=variable, y=value, fill=sources)) +
  geom_bar(stat="identity")+
  theme_bw(base_size = 30) + 
  theme(plot.title = element_text(size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title=element_text(face = 'bold'),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size=30),
        axis.ticks.x=element_blank()) + 
  top30_s_colors + 
  xlab("Participant") + 
  ylab("Dietary omega3 intake (g/day)") 
#scale_y_continuous(labels = percent_format(), expand = c(0, 0)) 
print(bar_cum)

ggsave(filename='./result/figures/distribution/fiber_subsets.png', plot=bar_cum, width = 15, height = 6, dpi = 600) 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Fig 1E. PCoA of species-level taxonomy (decorated with fruit fiber and CRP) ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
meta_pcoa <- mlvs_metadata_925[ , colnames(mlvs_metadata_925) %in% c('a_trn07_fo_dr','logcrp')] #a_omega3_fs_dr

# create quartiles for frtaf10v
quantile(meta_pcoa$a_trn07_fo_dr, na.rm=TRUE)
# 0%         25%         50%         75%        100% 
# 0.04428571  3.19000000  4.43000000  6.05571429 15.22000000 

#meta_pcoa <- meta_pcoa %>%
  #mutate(frtaf10vq = case_when(
    #0.789 < a_omega3_fs_dr & a_omega3_fs_dr <= 1.67 ~ 'Q1',
    #1.67 < a_omega3_fs_dr & a_omega3_fs_dr<= 2.21 ~ 'Q2',
    #2.21 < a_omega3_fs_dr & a_omega3_fs_dr<= 2.78 ~ 'Q3',
    #2.78 < a_omega3_fs_dr ~ 'Q4'))

meta_pcoa <- meta_pcoa %>%
  mutate(frtaf10vq = case_when(
    0.71 < a_trn07_fo_dr & a_trn07_fo_dr <= 1.30 ~ 'Q1',
    1.30 < a_trn07_fo_dr & a_trn07_fo_dr<= 1.80 ~ 'Q2',
    1.80 < a_trn07_fo_dr & a_trn07_fo_dr<= 2.41 ~ 'Q3',
    2.41 < a_trn07_fo_dr ~ 'Q4'))


table(meta_pcoa$frtaf10vq)
meta_pcoa$frtaf10vq<- factor(meta_pcoa$frtaf10vq, 
                             levels = c('Q1', 'Q2', 'Q3', 'Q4'))

# create quartiles for logcrp
quantile(meta_pcoa$logcrp, na.rm=TRUE)
# 0%        25%        50%        75%       100% 
# -1.6094379 -0.9162907 -0.1053605  0.6418539  2.9704145 

meta_pcoa <- meta_pcoa %>%
  mutate(logcrpq = case_when(
    -2.0 <= logcrp & logcrp <= -0.9162907 ~ 'Q1',
    -0.9162907 < logcrp & logcrp<= -0.1053605 ~ 'Q2',
    -0.1053605 < logcrp & logcrp<= 0.6418539 ~ 'Q3',
    0.6418539 < logcrp ~ 'Q4'))
table(meta_pcoa$logcrpq)

library(tidyverse)
meta_pcoa$logcrpq <- fct_explicit_na(meta_pcoa$logcrpq, na_level = "Missing")

rownames(meta_pcoa) <- rownames(mlvs_metadata_925)  



bugs_pcoa <- species_all
#make sure samples are in rows
bugs_pcoa <- capscale(bugs_pcoa ~ 1, distance = 'bray') #Question? old pcoa used asin sqrt transformed bug data?
variance.bugs <- head(eigenvals(bugs_pcoa)/sum(eigenvals(bugs_pcoa)))
x_variance.bugs <- as.integer(variance.bugs[1]*100)
y_variance.bugs <- as.integer(variance.bugs[2]*100)

# sample scores
ord.scores <- as.data.frame(scores(bugs_pcoa, display = "sites"))

# bug scores
ord.bug.scores <- as.data.frame(scores(bugs_pcoa, display = "species") )


# append metadata to scores
ord.scores.meta <- merge(ord.scores, meta_pcoa, by = 'row.names')

# set first column as rownames and remove it
rownames(ord.scores.meta) <- ord.scores.meta[ ,1]
ord.scores.meta[ ,1] <- NULL

# remove logCRP = NA
ord.scores.meta <- ord.scores.meta[complete.cases(ord.scores.meta), ]

# number of samples
n_samples <- nrow(ord.scores.meta)


library(viridis)
ordcolors <- scale_colour_gradientn(colours = viridis(10), limits=c(0,4))

pcoa.plot <- 
  ggplot( ord.scores.meta, aes(MDS1, MDS2) ) + 
  geom_point(size=3, alpha = 0.75, aes( color = a_trn07_fo_dr)) + #change name
  #shape = factor(logcrpq) or factor(a_omega3_fs_der)
  #scale_shape_manual(values=c(3, 16, 17, 15), name = "CRP quartiles")+
  labs(color = "Short Term Trans Fat (g/day)")+ #change name
  theme_classic() + 
  #theme_bw(base_size = 30) +
  theme(axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_text( size = 30, face = 'bold'), 
        legend.text = element_text(size=20), legend.title = element_text(size=30),
        legend.position = 'right', legend.direction = 'vertical', legend.box = 'vertical',
        plot.title = element_text( size = 30, face = 'bold')
  ) +
  guides( color = guide_legend(nrow=2, byrow = TRUE) ) +
  xlab(paste0("PCo1 (",x_variance.bugs,"%)")) + ylab(paste0("PCo2 (",y_variance.bugs,"%)"))

pcoa.plot <- pcoa.plot +   ordcolors + guides(colour="colourbar")

print(pcoa.plot)







# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### extra code ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#met <- read.csv('/Users/jorickbater/Downloads/mlvs_exposure_wenjie.csv', header=T, row.names = 1)
#met$id = rownames(met)

#z_metadata_all$choc_avg <- (z_metadata_all$dchocffq1+z_metadata_all$mchocffq1+z_metadata_all$dchocffq2+z_metadata_all$mchocffq2)/2
#z_metadata_all$choc_avg[is.na(z_metadata_all$choc_avg)] <- 0

#Blood Fatty acids
#fa <- readxl::read_xlsx("/Users/jorickbater/Desktop/Mingyang Stuff/fattyacids.xlsx")
#fa <- as.data.frame(fa)
#fa$ID <- as.numeric(fa$ID)
#fa <- fa[,3:47]
#dat <- merge(fa,tot,by="ID")

#gut <-read_tsv("metaphlan_taxonomic_profiles.tsv")
#gut <- t(gut)
#cols <- gut[1,]
#colnames(gut) <- cols
#gut = gut[-1,]
#gut <- gut[ , grepl( "s__" , colnames(gut) ) ]
#colnames(gut) <- gsub('.*s__', '', colnames(gut))

#gut <- as.data.frame(gut)
#gut <- setDT(gut, keep.rownames = TRUE)[]
#gut$rn <- gsub("([0-9]+)_.*", "\\1", gut$rn)

#gut <- sapply(gut[,1:544], as.numeric)
#gut <- as.data.frame(gut)
#gut$rn <- as.factor(gut$rn)
#gut <- aggregate(.~rn, gut, mean)
#gut <- rename(gut,ID=rn)

#rownames(gut) = gut$ID

#z_mlvs_exposure <- read.csv('/Users/jorickbater/Desktop/Mingyang Stuff/Data/mlvs_exposure_wenjie.csv', header=T, row.names = 1)
#z_mlvs_exposure$id = rownames(z_mlvs_exposure)
#z_mlvs_exposure$id <- as.factor(z_mlvs_exposure$id)
#z_mlvs_exposure <- merge(totom,z_mlvs_exposure,by="id")

#z_metadata_all$participant <- rownames(z_metadata_all)
#met <- z_metadata_all

#Merge ID
#ID <- read_sas("/Users/jorickbater/Desktop/Mingyang Stuff/Data/mlvs_bld_phase12.sas7bdat")
#ID = ID[,c("HarvardID","id")]

#tot <- merge(met,ID, by="id")
#tot$id <- NULL
#tot <- rename(tot,id=HarvardID)

#Gut Taxa
#tax_rpk_name <-   read.table( "./bugs_dna_929_unFilt.tsv",
#                              sep = '\t',    header = TRUE,    check.names = FALSE,    na.strings = c("", "NA"))
#dim(tax_rpk_name)
#tax_rpk_name<-tax_rpk_name %>%
#  separate(Sample, c("kingdom",       "phylum",        "class" ,        "order",         "family",        "genus" ,        "species" ,      "strain"), 
#           sep = '\\|', remove = TRUE)

# only keep species-level features
#tax_rpk_species <- subset(tax_rpk_name,!is.na(species) & is.na(strain))
#rownames(tax_rpk_species)<-tax_rpk_species$species
#tax_rpk_species<-tax_rpk_species[,-c(1:8)]

#all_species <- as.data.frame(tax_rpk_species)
#all_species <- all_species[ , !colnames(all_species) %in% c(grep("15074981", colnames(all_species), value = T))]
#dim(all_species)

#dim(all_species[apply(all_species, 1, function(bug) sum(bug >= 0.01) >= 0.1*ncol(all_species)), ])
#all_species_filt <- all_species[ apply(all_species, 1, function(bug) sum(bug >= 0.01) >= 0.1*ncol(all_species)), ]
#species_list <- rownames(all_species_filt)
#species_list.format <- tolower(gsub(pattern = "s__", "", species_list))

#z_idkey <- read_sas("/Users/jorickbater/Desktop/Mingyang Stuff/Data/idlinks.sas7bdat")
#z_idkey <- z_idkey[2:909,]
#rownames(z_idkey) <- z_idkey$id
#z_idkey <- z_idkey[order(z_idkey$id), ] 
#z_idkey <- rename(z_idkey,ID1=aliasid8digits)



#med_id<-subset(tot, select=id)
#ttax_rpk <- as.data.frame(t(tax_rpk_species))
#ttax_rpk$id<-rownames(ttax_rpk)
#ttax_rpk$id <- gsub("([0-9]+)_.*", "\\1", ttax_rpk$id)

#ttax_rpk<-inner_join(med_id, ttax_rpk, by="id")
#rownames(ttax_rpk)<-ttax_rpk$id
#ttax_rpk <- subset(ttax_rpk, select = -id)
#View(ttax_rpk)

#dat <- tot
#nam = as.vector(dat$ID)
#gut <- filter(gut, ID %in% nam)

#gut$ID <- NULL
#gut[,2:544] <- gut[,2:544]/100

#rownames(dat) = dat$ID
#dat$ID = as.factor(dat$ID)

#Make Phyloseq Object
#OTU = t(gut)
#OTU = otu_table(OTU, taxa_are_rows = TRUE)
#MET = sample_data(dat, errorIfNULL=FALSE)
#phy <- phyloseq(OTU, MET)

#Filter
#library(metagMisc)
#phy <- phyloseq_filter_prevalence(phy, prev.trh = 0.10, abund.trh = 0.01, threshold_condition = "AND") 
#gut = t(phy@otu_table@.Data)
#dim(gut)

#Alpha Diversity
#plot_richness(phy, x="P3", measures = c("Shannon", "InvSimpson")) + geom_boxplot()
#plot_richness(phy, x="P10", measures = c("Shannon", "InvSimpson")) + geom_boxplot()
#plot_richness(phy, x="P13", measures = c("Shannon", "InvSimpson")) + geom_boxplot()
#plot_richness(phy, x="P14", measures = c("Shannon", "InvSimpson")) + geom_boxplot()
#plot_richness(phy, x="crp_plasma1", measures = c("Shannon", "InvSimpson")) + geom_boxplot()
#plot_richness(phy, x="omega_3_epa_ffq1", measures = c("Shannon", "InvSimpson")) + geom_boxplot()

#Beta Diversity
#source("miseqR.R")
#set.seed(1)
#dat_scale <- phy %>%
#  scale_reads(round = "round") 
#dat_bray <- phyloseq::distance(dat_scale, method = "bray")
#ordination = ordinate(phy, method="PCoA", distance=dat_bray)
#plot_ordination(dat, ordination, color="SampleType") + theme(aspect.ratio=1)

#PERMANOVA
#library(microbiome)
#library(vegan)
#rel <- microbiome::transform(phy, "compositional")

#Beta Diversity

#p <- plot_landscape(rel, method = "NMDS", distance = "bray", col = "crp_avg", size = 3)
#plot(p)

#otu <- abundances(phy)
#meta <- meta(phy)
#permanova <- adonis(t(otu) ~ P3 + P10 + P13 + P14 + crp_avg,
#data = dat, permutations=99, method = "bray")
#permanova <- adonis(t(otu) ~ crp_avg,
 #                   data = meta, permutations=999, method = "bray")

#permanova

#coef <- coefficients(permanova)["crp_avg",]
#top.coef <- coef[rev(order(abs(coef)))[1:20]]
#par(mar = c(3, 14, 2, 1))
#barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")

#Run Maaslin2
#library(Maaslin2)
#fit_data <- Maaslin2(
  #gut, dat, 'demo_output', transform = "none", heatmap_first_n = 100,
  #max_significance = 0.25,
  #fixed_effects = c('cod_liv_oil_avg','linseed_avg','flax_avg','omega_3_epa_avg', 'walnuts_avg', 'avocado_avg', 'leafvegavg', 'fishavg', 'agemlvs','calor_avg','abx_avg'),
  #random_effects = 'ID1',
  #normalization = 'NONE',
  #standardize = FALSE)


#FFQ Data
#ffq1 <- read_sas("/Users/jorickbater/Desktop/Mingyang Stuff/Data/mlvsffq1_excl_nts.sas7bdat")
#om1 <- ffq1[,c("id","pfa183n3c07_fs_ffq1", "pf205n3c07_fs_ffq1", "pf225n3c07_fs_ffq1", "pf226n3c07_fs_ffq1")]
#ffq2 <- read_sas("/Users/jorickbater/Desktop/Mingyang Stuff/Data/mlvsffq2_excl_nts.sas7bdat")
#om2 <- ffq2[,c("id","pfa183n3c07_fs_ffq2", "pf205n3c07_fs_ffq2", "pf225n3c07_fs_ffq2", "pf226n3c07_fs_ffq2")]
#dr1 <- read_sas("/Users/jorickbater/Desktop/Mingyang Stuff/Data/dr_nts_wk1_mean_e_adjusted.sas7bdat")
#dom1 <- dr1[,c("id","a_omega3_fs_dr_w1avg", "a_p22_5_fs_dr_w1avg")]
#dr2 <- read_sas("/Users/jorickbater/Desktop/Mingyang Stuff/Data/dr_nts_wk2_mean_e_adjusted.sas7bdat")
#dom2 <- dr2[,c("id","a_omega3_fs_dr_w2avg", "a_p22_5_fs_dr_w2avg")]

#dat_list = list(om1,om2,dom1,dom2)

#my_merge <- function(df1, df2){                               
#  merge(df1, df2, by = "id")
#}

#totom <- Reduce(my_merge,dat_list)
#totom$pfa183n3c07_fs_avg <- rowMeans(totom[,c('pfa183n3c07_fs_ffq1', 'pfa183n3c07_fs_ffq2')], na.rm=TRUE)
#totom$pf205n3c07_fs_avg <- rowMeans(totom[,c('pf205n3c07_fs_ffq1', 'pf205n3c07_fs_ffq2')], na.rm=TRUE)
#totom$pf225n3c07_fs_avg <- rowMeans(totom[,c('pf225n3c07_fs_ffq1', 'pf225n3c07_fs_ffq2')], na.rm=TRUE)
#totom$pf226n3c07_fs_avg <- rowMeans(totom[,c('pf226n3c07_fs_ffq1', 'pf226n3c07_fs_ffq2')], na.rm=TRUE)

