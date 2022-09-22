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
source("noGit/Rstart.R")





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
########################## (1) metadata Exposure Data ###############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

###################### read in channing metadata ###############################


z_mlvs_exposure_new <- read.table("noGit/mlvs_exposure_for_jorick.txt", 
                                  header=TRUE, sep='\t', check.names=TRUE, quote ="")
#z_mlvs_exposure_new <- merge(totom,z_mlvs_exposure_new,by="id")

z_idkey <- read.csv('noGit/idkey.csv', header=TRUE, sep = ",", stringsAsFactors = FALSE)
rownames(z_idkey) <- z_idkey$id
z_idkey <- z_idkey[order(z_idkey$id), ] 
z_mlvs_exposure <- reassignKey(z_mlvs_exposure_new)
z_mlvs_exposure$id <- as.integer(rownames(z_mlvs_exposure))


########+ Plasma Data  from Mingyang#########
setwd("..")
z_plasma <- read.table("noGit/plasma_fatty_acids_sjbae.txt",
                       header=TRUE, sep='\t', check.names =TRUE, quote ="")
z_plasma$id <- as.integer(substr(z_plasma$ID, 1, 6)) 
#get first 6 number of ID1 which stands for id in idkey file

z_mlvs_exposure_plasma<- left_join(z_mlvs_exposure, z_plasma, by ="id")
write.table(z_mlvs_exposure_plasma, file = "noGit/matching_exposure_plasma.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)







#gsub search and replace in files -> ~~_w1avg colnames changed to ~~_w1, etc.
colnames(z_mlvs_exposure_plasma) <- gsub("w1avg","w1", colnames(z_mlvs_exposure_plasma)) 
colnames(z_mlvs_exposure_plasma) <- gsub("w2avg","w2", colnames(z_mlvs_exposure_plasma))
colnames(z_mlvs_exposure_plasma) <- gsub("plasma1","w1", colnames(z_mlvs_exposure_plasma))
colnames(z_mlvs_exposure_plasma) <- gsub("plasma2","w2", colnames(z_mlvs_exposure_plasma))

z_mlvs_agebmi <- z_mlvs_exposure_plasma [ , colnames(z_mlvs_exposure_plasma) %in% c('agemlvs', 'bmi12','act10','id')] 
#id is actually rownames

#cohort FFQ, ignore for now
keep_metadata_cum <- c('ala10v','epa10v','dha10v','dpa10v','calor10v', 'trans10v', 'omega610v', 'omega310v', 'omega3_noala10v')
z_metadata_cum <- z_mlvs_exposure_plasma[, keep_metadata_cum]

keep_metadata_avg <- c(
  #ffq questionnaire
  'acid_avg', 'alc_avg', 'abx_avg', 'prep_avg',
  'selfdiet_avg', 'probx_avg', 'bristolcat_avg', 'bristol_avg', 'oral_avg', 'yogurt_avg',
  #ffq energy
  'calor_avg',
  #Omega nutrients
  'ala_avg', 'epa_avg', 'dha_avg', 'dpa_avg', 'omega3_avg', 'omega3_noala_avg', 'omega6_avg','trans_avg'
)

z_metadata_avg <- z_mlvs_exposure_plasma[, keep_metadata_avg]

keep_metadata_w1 <- c(
  #questionnaire
  'acid_w1', 'alc_w1', 'abx_w1', 'prep_w1',
  'selfdiet_w1', 'probx_w1', 'bristolcat_w1', 'bristolcat5', 'bristolcat6', 'oral_w1', 'yogurt_w1',
  'bristol_w1', 'bristol5', 'bristol6',
  #ddr1
  'calor_fs_dr_w1', 
  #ddr1 nutrients
  #these are blood????
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

z_metadata_w1 <- z_mlvs_exposure_plasma[, keep_metadata_w1]
z_metadata_w2 <- z_mlvs_exposure_plasma[, keep_metadata_w2]

keep_plasmafattyacid <- c('ala_pfa', 'epa_pfa', 'dpa_pfa', 'dha_pfa')


z_metadata_pfa <- z_mlvs_exposure_plasma[, keep_plasmafattyacid]

z_metadata_all <- cbind(z_mlvs_agebmi, z_metadata_cum, z_metadata_avg, z_metadata_w1, z_metadata_w2, z_metadata_pfa)
#z_metadata_all <- cbind(z_mlvs_agebmi, z_metadata_avg, z_metadata_w1, z_metadata_w2)


##### + subset metadata to prep for maaslin.pcl files ######

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
############################## (3) species data ###############################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# bring in legacy file - just for checking, I actually used species data in Channing
# library(dplyr)
all_species <- read_tsv("noGit/metaphlan_taxonomic_profiles.tsv")
names(all_species)[names(all_species) == '# taxonomy'] <- 'Sample'

all_species <- all_species %>% 
  rename_at(.vars = vars(ends_with("_dna_taxonomic_profile")),
            .funs = funs(sub("[_]dna_taxonomic_profile", "", .)))

all_species <-all_species %>%
  separate(Sample, c("kingdom","phylum","class" ,"order","family","genus" ,"species" ), 
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



###### TAX and OTU table
TAX <- as.matrix(all_species1[,1:7])

OTU <- all_species1[,8:ncol(all_species1)]
# OTU maybe later

######

all_species1<-all_species1[,-c(1:7)]
#all_species1<-all_species1[,-c(1:8)]

all_species1 <- as.data.frame(all_species1)

dim(all_species1)
# [1] 468 929
# sunjeong: [1] 543 929

all_sample <- colnames(all_species1)
#View(all_sample)

# remove guy with colectomy 
all_species1 <- all_species1[ , !colnames(all_species1) %in% c(grep("15074981", colnames(all_species1), value = T))]
dim(all_species1)
# [1] 468 925
# sunjeong: [1] 543 925
#colSums(all_species)

# generate list of species after abundance/prevalence filtering
# since this file is scaled to 100, i am keeping bugs with >10% prevalence and 0.01 abundance (i.e. .0001*100 when relab is scaled 0-1 and not 0-100)
dim(all_species1[apply(all_species1, 1, function(bug) sum(bug >= 0.01) >= 0.1*ncol(all_species1)), ])
# [1] 139 925 = the 139 species after QC included in initial starr manuscripts
# sunjeong: [1] 156 925
all_species_filt <- all_species1[ apply(all_species1, 1, function(bug) sum(bug >= 0.01) >= 0.1*ncol(all_species1)), ]





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# important! bug tables from channing need to be converted to 0-1 (not 0-11), which will wreck our downstream stats

all_species_filt   <- all_species_filt  / 100 

# let's remove all the stool suffixes from the sas names to facilitate changing the name back to species names 
#colnames(z_species.s1) <- gsub("*_s\\d","", colnames(z_species.s1))

# the species i want to keep are designated by filtering above (species_list), gonna lower case it and get rid of leading "s__"
#_rownames(all_species_filt) <- tolower(gsub(pattern = "s__", "", rownames(all_species_filt)))
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
z_metadata_all$participant <- z_metadata_all$id

# merge the species and metadata tables 
z_merged.all <- merge(x = z_species.all, y = z_metadata_all, by = 'participant', all.x = TRUE)

# rename the columns so they're unique participant_s#
rownames(z_merged.all) <- paste(z_merged.all$participant, z_merged.all$stoolnum, sep='_s') #naming convention

#checkit <- colnames(z_merged.all)
#View(checkit)


# keep only related metadata to the time course we are interested in 
# reorder so that metadata on top, bacterial data on the bottom 
z_ffq.all_species.all <- z_merged.all %>% select(participant, stoolnum, colnames(z_mlvs_agebmi), 
                                                 colnames(z_metadata_avg), everything(), 
                                                 -contains('w1'), -contains('w2'))
z_w1.all_species.all <- z_merged.all %>% select(participant, stoolnum, colnames(z_mlvs_agebmi), 
                                                colnames(z_metadata_w1), contains ('w1'), everything(), 
                                                -contains('avg'), -contains('w2'))
z_w2.all_species.all <- z_merged.all %>% select(participant, stoolnum, colnames(z_mlvs_agebmi), 
                                                colnames(z_metadata_w2), contains ('w2'), everything(), 
                                                -contains('avg'), -contains('w1'))

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






################ [2] metadata set #########


###### (1) cohort FFQ cumulative average ######
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

#### (2)  MLVS FFQ#######
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

###### (3) MLVS DDR #####
#'a_omega3_fs', 'a_omega3_noala_fs', 'a_ala_fs_dr_w1', 'a_f205_fs_dr_w1', 'a_f226_fs_dr_w1', 'a_p22_5_fs_dr_w1', 'a_trn07_fo_dr_w1'
ddr_keep <- c('sample','participant', 'agemlvs', 'calor_fs_dr', 'abx')
ddr_keep_nocal <- c('sample','participant', 'agemlvs', 'abx')

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












#### (4) MLVS biomarker ######
biomarker_keep <- c('sample','participant', 'agemlvs', 'abx', 'bmi12', 'calor_fs_dr')
biomarker_logcrp  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logcrp')] 
biomarker_logcrp_complete <- biomarker_logcrp[!is.na(biomarker_logcrp$logcrp), ]
biomarker_loghdl  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'loghdl')] 
biomarker_loghdl_complete <- biomarker_loghdl[!is.na(biomarker_loghdl$loghdl), ]
biomarker_logtg  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtg')] 
biomarker_logtg_complete <- biomarker_logtg[!is.na(biomarker_logtg$logtg), ]
## just logcrp vs complete just the same for maaslin




#### biomarker and omega3 ####

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


write.table(mlvs_metadata_925, file = "noGit/mlvs_metadata_925.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#############  [3] pcl writing  ############
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

for(i in 1:ncol(species_all)) { 
  species_all[ , i] <- as.numeric(as.character(species_all[, i])) 
}

pcl_file <- c('species_all', 
              'ala_ffqcum','epa_ffqcum','dha_ffqcum', 'dpa_ffqcum', 'trans_ffqcum', 'omega6_ffqcum', 'omega3_ffqcum', 'omega3_noala_ffqcum',
              'ala_avg_fs_ffq', 'epa_avg_fs_ffq', 'dha_avg_fs_ffq', 'dpa_avg_fs_ffq', 'trans_avg_fs_ffq', 'omega6_avg_fs_ffq', 'omega3_avg_fs_ffq', 'omega3_noala_avg_fs_ffq',
              'a_omega3_fs_ddr', 'a_omega3_fs_ddr_complete', #'a_omega3_fs_ddr_1.5IQR', 'a_omega3_fs_ddr_3IQR', 
              'a_omega3_noala_fs_ddr', 'a_ala_fs_ddr', 'a_trans_fs_ddr', 'a_epa_fs_ddr', 'a_dha_fs_ddr', 'a_dpa_fs_ddr',
              'omega3_tfat_ddr', 'trans_tfat_ddr', 'epa_tfat_ddr', 'dha_tfat_ddr', 'dpa_tfat_ddr',
              'biomarker_logcrp', 'biomarker_logcrp_complete', #'a_omega3_fs_ddr_logcrp', 
              'biomarker_loghdl', 'biomarker_loghdl_complete', #'a_omega3_fs_ddr_loghdl', 
              'biomarker_logtg', 'biomarker_logtg_complete', #'a_omega3_fs_ddr_logtg',
              #'a_omega3_noala_fs_ddr_logcrp', 'a_omega3_noala_fs_ddr_loghdl', 'a_omega3_noala_fs_ddr_logtg',
              # 'a_trans_fs_ddr_logcrp', 'a_trans_fs_ddr_loghdl', 'a_trans_fs_ddr_logtg',
              # 'a_ala_fs_ddr_logcrp', 'a_ala_fs_ddr_loghdl', 'a_ala_fs_ddr_logtg',
              # 'a_epa_fs_ddr_logcrp', 'a_epa_fs_ddr_loghdl', 'a_epa_fs_ddr_logtg',
              # 'a_dha_fs_ddr_logcrp', 'a_dha_fs_ddr_loghdl', 'a_dha_fs_ddr_logtg',
              # 'a_dpa_fs_ddr_logcrp', 'a_dpa_fs_ddr_loghdl', 'a_dpa_fs_ddr_logtg',
              'z_metadata_all','mlvs_metadata_925')   

# make pcl files 
for (a in 1:length(pcl_file)){
  tempdf <- NULL
  # transpose
  tempdf <- as.data.frame ( t (get( pcl_file[a] ) ) )
  # steps to add the word sample to cell [1,1] so maaslin won't throw an error 
  tempcol <- colnames(tempdf)
  tempdf$sample = row.names(tempdf)
  addsamp = c('sample', tempcol)
  foo = tempdf [ , addsamp]
  # write table
  write.table(foo, paste('./noGit/maaslin2/pcl/', pcl_file[a], '.pcl',sep=''),  sep = '\t', quote = F, eol = '\n',  row.names=F)
} 

## Done preparing for analyses
# now you can start Analysis.R (even with empty enviornment)




