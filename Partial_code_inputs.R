#####Partial code from Jorick including input data#####


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
########################## prepare dataset ####################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
library(haven)
library(data.table)
library(dplyr)
library(tidyverse)

setwd("/Users/jorickbater/Downloads")

### Pre-processing steps 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
########################## (1) metadata ###############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

###################### read in channing metadata ###############################

#Exposure Data
source("/Users/jorickbater/Desktop/Mingyang Stuff/Data/Rstart.R")
z_mlvs_exposure_new <- read.table('./mlvs_exposure_for_jorick.txt', header=TRUE, sep='\t', check.names=TRUE, quote ="")
#z_mlvs_exposure_new <- merge(totom,z_mlvs_exposure_new,by="id")

z_idkey <- read.csv('/Users/jorickbater/Desktop/Mingyang Stuff/Data/idkey.csv', header=TRUE, sep = ",", stringsAsFactors = FALSE)
rownames(z_idkey) <- z_idkey$id
z_idkey <- z_idkey[order(z_idkey$id), ] 

z_mlvs_exposure <- reassignKey(z_mlvs_exposure_new)

colnames(z_mlvs_exposure) <- gsub("w1avg","w1", colnames(z_mlvs_exposure))
colnames(z_mlvs_exposure) <- gsub("w2avg","w2", colnames(z_mlvs_exposure))

colnames(z_mlvs_exposure) <- gsub("plasma1","w1", colnames(z_mlvs_exposure))
colnames(z_mlvs_exposure) <- gsub("plasma2","w2", colnames(z_mlvs_exposure))

z_mlvs_agebmi <- z_mlvs_exposure [ , colnames(z_mlvs_exposure) %in% c('agemlvs', 'bmi12','act10','id')] 

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