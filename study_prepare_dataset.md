# [1] prepare dataset
~~~~~~~~~~~r~~~~~~~~~~~~~
library(haven)
library(data.table)
library(dplyr)
library(tidyverse)
~~~~~~~~~~~~~~~~~~~~~~~~
libary used
source("noGit/Rstart.R")

## (1) metadata
### read in channing metadata 

```diff
- what is channing metadata??
```

1) Exposure Data
~~~~~~~~~~~r~~~~~~~~~~~~~
z_mlvs_exposure_new <- read.table("noGit/mlvs_exposure_for_jorick.txt", 
                                  header=TRUE, sep='\t', check.names=TRUE, quote ="")
~~~~~~~~~~~~~~~~~~~~~~~~
these data seems to have all the information about the subjects.
including height, age, smoking, bmi, etc.

2) ID
~~~~~~~~~~~r~~~~~~~~~~~~~
z_idkey <- read.csv('noGit/idkey.csv', header=TRUE, sep = ",", stringsAsFactors = FALSE)
rownames(z_idkey) <- z_idkey$id
z_idkey <- z_idkey[order(z_idkey$id), ] 
z_mlvs_exposure <- reassignKey(z_mlvs_exposure_new)
~~~~~~~~~~~~~~~~~~~~~~~~
```diff
- suppose reassigning key for blinding and/or randomization?
```

3) Change names
gsub search and replace in files -> ~_w1avg colnames changed to ~_w1, etc.
```diff
- suppose w1avg is week 1 average
```
~~~~~~~~~~r~~~~~~~~~~~~~~
colnames(z_mlvs_exposure) <- gsub("w1avg","w1", colnames(z_mlvs_exposure)) 
colnames(z_mlvs_exposure) <- gsub("w2avg","w2", colnames(z_mlvs_exposure))
colnames(z_mlvs_exposure) <- gsub("plasma1","w1", colnames(z_mlvs_exposure))
colnames(z_mlvs_exposure) <- gsub("plasma2","w2", colnames(z_mlvs_exposure))
~~~~~~~~~~~~~~~~~~~~~~~~


4) get age bmi -> basic informations
~~~~~~~~~~r~~~~~~~~~~~~~~
z_mlvs_agebmi <- z_mlvs_exposure [ , colnames(z_mlvs_exposure) %in% c('agemlvs', 'bmi12','act10','id')] 
#id is actually rownames
~~~~~~~~~~~~~~~~~~~~~~~~
5) extracted omega3 related info -> interventions
~~~~~~~~~~r~~~~~~~~~~~~~~
keep_metadata_cum <- c('ala10v','epa10v','dha10v','dpa10v','calor10v', 'trans10v', 'omega610v', 'omega310v', 'omega3_noala10v')
z_metadata_cum <- z_mlvs_exposure[, keep_metadata_cum]
~~~~~~~~~~~~~~~~~~~~~~~~
```diff
- what does 10v stands for??
```
6) extract average
~~~~~~~~~~r~~~~~~~~~~~~~~
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
~~~~~~~~~~~~~~~~~~~~~~~~
```diff
- Who's average?? total or single participant.
- probably single participant
- then what type of average?
```
7) week 1, week 2
~~~~~~~~~~~r~~~~~~~~~~~~~
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
~~~~~~~~~~~~~~~~~~~~~~~~
9) add up all the wanted data
z_metadata_all <- cbind(z_mlvs_agebmi, z_metadata_cum, z_metadata_avg, z_metadata_w1, z_metadata_w2)
#z_metadata_all <- cbind(z_mlvs_agebmi, z_metadata_avg, z_metadata_w1, z_metadata_w2)

## subset metadata to prep for maaslin.pcl files 
```diff
- When do we use this subset?
```
~~~~~~~~~~~r~~~~~~~~~~~~~
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
~~~~~~~~~~~~~~~~~~~~~~~~

# (2) species data
```diff
- What is channing.... 
```
~~~~~~~~~~~~r~~~~~~~~~~~~
# bring in legacy file - just for checking, I actually used species data in Channing
all_species <- read_tsv("noGit/metaphlan_taxonomic_profiles.tsv")
names(all_species)[names(all_species) == '# taxonomy'] <- 'Sample'
~~~~~~~~~~~~~~~~~~~~~~~~

change the names. erasing "_dna_taxonomic_profile"
~~~~~~~~~~~r~~~~~~~~~~~~~
all_species <- all_species %>% 
  rename_at(.vars = vars(ends_with("_dna_taxonomic_profile")),
            .funs = funs(sub("[_]dna_taxonomic_profile", "", .)))
~~~~~~~~~~~~~~~~~~~~~~~~
separate taxanomy
~~~~~~~~~~~~~~~~~~~~~~~~
all_species <-all_species %>%
  separate(Sample, c("kingdom","phylum","class" ,"order","family","genus" ,"species" ), 
           sep = '\\|', remove = TRUE)
~~~~~~~~~~~~~~~~~~~~~~~~

leave only species-level features
~~~~~~~~~~~~r~~~~~~~~~~~~
all_species1 <- subset(all_species,!is.na(species))
all_species1 <- as.data.frame(all_species1) #weird formatting change to maintain rows
rownames(all_species1)<-all_species1$species
all_species1<-all_species1[,-c(1:7)]
all_species <- as.data.frame(all_species1)
~~~~~~~~~~~~~~~~~~~~~~~~


# take out un-relevant features

### patient outlier

how many species?
~~~~~~~~~~~r~~~~~~~~~~~~~
dim(all_species)
all_sample <- colnames(all_species)
~~~~~~~~~~~~~~~~~~~~~~~~
sunjeong: [1] 543 929

remove guy with colectomy 
~~~~~~~~~~~r~~~~~~~~~~~~~
all_species <- all_species[ , !colnames(all_species) %in% c(grep("15074981", colnames(all_species), value = T))]
dim(all_species)
~~~~~~~~~~~~~~~~~~~~~~~~
sunjeong: [1] 543 925 - this means we had 4 samples from colectomy patient

### abundance/prevalence filtering
since this file is scaled to 100, i am keeping bugs with >10% prevalence and 0.01 abundance (i.e. .0001*100 when relab is scaled 0-1 and not 0-100)
0.0001 이상으로 존재하는 샘플의 수가 10% 이상 이라는 뜻
~~~~~~~~~~~r~~~~~~~~~~~~~
all_species_filt <- all_species[ apply(all_species, 1, function(bug) sum(bug >= 0.01) >= 0.1*ncol(all_species)), ]
dim(all_species_filt)
~~~~~~~~~~~~~~~~~~~~~~~~~
sunjeong: [1] 156 925
after QC included in initial starr manuscripts


important! bug tables from channing need to be converted to 0-1 (not 0-11), which will wreck our downstream stats
```diff
- What is "(not 0-11)"? Isn't it should be something like "(not 0-100%)"??
```
~~~~~~~~~~~r~~~~~~~~~~~~~
all_species_filt   <- all_species_filt  / 100 
~~~~~~~~~~~~~~~~~~~~~~~~

# remove all the stool suffixes from the sas names
1) species 
the species i want to keep are designated by filtering above (species_list), gonna lower case it and get rid of leading "s__"
```diff
- What is "(species_list)"? Isn't it should be something like "(all_species_filt)"??
```
~~~~~~~~~~~r~~~~~~~~~~~~~
rownames(all_species_filt) <- tolower(gsub(pattern = "s__", "", rownames(all_species_filt)))
all_species_filt <- t(all_species_filt)
all_species_filt <- as.data.frame(all_species_filt)
~~~~~~~~~~~~~~~~~~~~~~~~

2) ID suffixes
~~~~~~~~~~~r~~~~~~~~~~~~~
all_species_filt$id <- substr(rownames(all_species_filt), start = 1, stop = 6)
all_species_filt$st <- substr(rownames(all_species_filt), start = 10, stop = 13)
all_species_filt$st <- str_replace_all(all_species_filt$st,  c("SF05" = "s1","SF06" = "s2","SF07" = "s3","SF08" = "s4"))
rownames(all_species_filt) <- paste0(all_species_filt$id,"_", all_species_filt$st)
~~~~~~~~~~~~~~~~~~~~~~~~

3) change name id st
~~~~~~~~~~~r~~~~~~~~~~~~~
z_species.all = all_species_filt
names(z_species.all)[names(z_species.all) == 'id'] <- 'participant'
names(z_species.all)[names(z_species.all) == 'st'] <- 'stoolnum'
z_species.all$stoolnum <- str_replace_all(z_species.all$stoolnum,  c("s" = ""))
species_all <- all_species_filt[,!colnames(all_species_filt) %in% c('id','st')]
~~~~~~~~~~~~~~~~~~~~~~~~






# append metadata onto short- and long-term .pcl files

1) metadata
metadata_all doesn't have a participants variable 
~~~~~~~~~~~r~~~~~~~~~~~~~
z_metadata_all$participant <- rownames(z_metadata_all)
~~~~~~~~~~~~~~~~~~~~~~~~

merge the species and metadata tables 
~~~~~~~~~~~r~~~~~~~~~~~~~
z_merged.all <- merge(x = z_species.all, y = z_metadata_all, by = 'participant', all.x = TRUE)
~~~~~~~~~~~~~~~~~~~~~~~~

rename the columns so they're unique participant_s#
~~~~~~~~~~~r~~~~~~~~~~~~~
rownames(z_merged.all) <- paste(z_merged.all$participant, z_merged.all$stoolnum, sep='_s') #naming convention
~~~~~~~~~~~~~~~~~~~~~~~~


keep only related metadata to the time course we are interested in 
reorder so that metadata on top, bacterial data on the bottom 
~~~~~~~~~~~~r~~~~~~~~~~~~
z_ffq.all_species.all <- z_merged.all %>% select(participant, stoolnum, colnames(z_mlvs_agebmi), 
                                                 colnames(z_metadata_avg), everything(), 
                                                 -contains('w1'), -contains('w2'))
z_w1.all_species.all <- z_merged.all %>% select(participant, stoolnum, colnames(z_mlvs_agebmi), 
                                                colnames(z_metadata_w1), contains ('w1'), everything(), 
                                                -contains('avg'), -contains('w2'))
z_w2.all_species.all <- z_merged.all %>% select(participant, stoolnum, colnames(z_mlvs_agebmi), 
                                                colnames(z_metadata_w2), contains ('w2'), everything(), 
                                                -contains('avg'), -contains('w1'))
~~~~~~~~~~~~~~~~~~~~~~~~                                                

2) subset stool
~~~~~~~~~~~~r~~~~~~~~~~~~  
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

~~~~~~~~~~~~~~~~~~~~~~~~

3) combime by week
finally, you've got one dataset with w1 stool and w1 metadata and another with w2 stool with w2 metadata 
now we can combine them back into one, large flat table to do the short-term analysis but one barrier is that the metadata in each 
has a _w1 or _w2 suffix. to fix this, rename the metadata variables to remove this so you can rbind 
~~~~~~~~~~~~r~~~~~~~~~~~~
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

z_w12.all_species.all <- rbind (z_s1.all_species.all, z_s2.all_species.all, z_s3.all_species.all, z_s4.all_species.all)
~~~~~~~~~~~~~~~~~~~~~~~~  
```diff
- Why are we binding different time points to same names?
```
BEAUTIFUL! now you have full metadata :: species tables for both long- and short-term

1. z_ffq.all_species.all = long-term :: species
2. z_w12.all_species.all = short-term :: species 






