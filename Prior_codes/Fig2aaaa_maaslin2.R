
##############################################################################################################################################
# 1) Purposes: conduct MaAsLin regression to calculate the associations of the measures of physical activity and adiposity with abundances of per microbial species
# 2) Study design:  Cross-sectional design
# 3) Endpoints: primary: Long-term body weight change from age 21 to the time of stool collection
secondary: BMI and fat mass% at stool collection, body weight change between the two collections of stool sample, and biomarkers of HbA1c and CRP
# 4) Exposures: primary: cumulative average of long-term physical activity level
secondary: short-term physical activity level measured by accelerometer, long- and short-term physical activity by intensity
# 5) Covariates: age, diet quality, total energy intake, smoking status, antibiotic use, probiotic use, stool type
# 6) Follow-up: 2012 (Men's Lifestyle Validation Study)  
#############################################################################################################################################

rm(list=ls())
getwd()
setwd("/noGit")
options(max.print=1000000)

library(Maaslin2)
library(ggplot2)
library(grid)
library(ggpubr)
library(tidyverse)
library(gridExtra)

#### Import data ####
# read in metadata
totmeta_stn <- read.table('./data/totmeta_stn.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")

# read in taxonomy data
species <-read.table(file = './data/species_filt.csv',
                     sep = ',',  row.names = 1,  header = TRUE,   check.names = FALSE)

# maaslin regressions
maaslin_diet_ddr<-function(exposure){
  dir = paste('./maaslin_results/', exposure, '.pcl',sep='')
  Maaslin2(species, 
           totmeta_stn, 
           dir, ##### input
           transform="AST",
           min_abundance=0.0001,
           min_prevalence=0.1,
           random_effects="ID1", 
           fixed_effects = c(exposure, "agemlvs","calor_fs_dr_ddr", "probx_ddr", "abx_ddr", "bristol_ddr", "smoke12" ))
}


maaslin_diet_ddr(exposure = 'wtchg')
maaslin_diet_ddr(exposure = 'logomega310v')

for (e in exposure_long){
  maaslin_diet_ddr(exposure = e)
} #doing

exposure_short <- c( 'ala_ddr','epa_ddr' ,'dha_ddr','dpa_ddr' , 'omega3_noala_ddr' , 'omega3_ddr',  'a_fat_fs_dr_ddr') 
for (e in exposure_short){
  maaslin_diet_ddr(exposure = e)
}

for (e in exposure_1yr){
  maaslin_diet_ddr(exposure = e)
}

for (e in exposure_avg){
  maaslin_diet_ddr(exposure = e)
}

for (e in outcome_plasma){
  maaslin_diet_ddr(exposure = e)
} #done

outcome_pfa <- c('ala_pfa', 'epa_pfa', 'dpa_pfa', 'dha_pfa',  'AA_pfa',  'omega3_pfa','omega3_noala_pfa','omega6_pfa')
for (e in outcome_pfa){
  maaslin_diet_ddr(exposure = e)
}

fatty_acid <- c('sat_pfa',
                'omega6_pfa','sat_pfa', 'monounsat_pfa', 'trans_pfa')
for (e in fatty_acid){
  maaslin_diet_ddr(exposure = e)
}

obesity <- c( 'wt_paq', 'height_paq1', 'bmi_paq1', 'waist_paq1', 'hip_paq1', 'whr_paq1', 'wtchg')
for (e in obesity){
  maaslin_diet_ddr(exposure = e)
}

maaslin_diet_ddr(exposure = 'wtchg')
### log it

logall <- c(exposure_long, exposure_short, exposure_1yr, exposure_avg, outcome_plasma, outcome_pfa)

logall <- paste("log", logall, sep = "")


for (e in logall){
  maaslin_diet_ddr(exposure = e)
}



##### maaslin of pfa adjusted by ddr #####
adjust <- c('ala', 'epa', 'dpa', 'dha','omega3','omega3_noala','omega6')
adjust <- c('omega3',
                'omega6','trans_pfa')

maaslin_adjust<-function(adjust){
  dir = paste('./maaslin_results/adjusted_', adjust, '.pcl',sep='')
  Maaslin2(species, 
           totmeta_stn, 
           dir, ##### input
           transform="AST",
           min_abundance=0.0001,
           min_prevalence=0.1,
           random_effects="ID1", 
           fixed_effects = c(paste(adjust, "_pfa", sep = ""),paste(adjust, "_ddr", sep = ""), "agemlvs","calor_fs_dr_ddr", "probx_ddr", "abx_ddr", "bristol_ddr", "smoke12" ))
}

for (a in adjust){
  maaslin_adjust(adjust = a)
}

##### maaslin of interaction term #####
adjust <- c('ala', 'epa', 'dpa', 'dha','omega3','omega3_noala','omega6')

maaslin_adjust<-function(adjust){
  dir = paste('./maaslin_results/adjusted_', adjust, '.pcl',sep='')
  Maaslin2(species, 
           totmeta_stn, 
           dir, ##### input
           transform="AST",
           min_abundance=0.0001,
           min_prevalence=0.1,
           random_effects="ID1", 
           fixed_effects = c(paste(adjust, "_pfa", sep = ""),paste(adjust, "_ddr", sep = ""), "agemlvs","calor_fs_dr_ddr", "probx_ddr", "abx_ddr", "bristol_ddr", "smoke12" ))
}

for (a in adjust){
  maaslin_adjust(adjust = a)
}


##### maaslin of interaction term intake/pfa #####
dir = paste('./maaslin_results/adjusted_ala_to_omega3.pcl',sep='')
Maaslin2(species, 
         totmeta_stn, 
         dir, ##### input
         transform="AST",
         min_abundance=0.0001,
         min_prevalence=0.1,
         random_effects="ID1", 
         fixed_effects = c('ala_ddr', 'omega3_noala_pfa', "agemlvs","calor_fs_dr_ddr", "probx_ddr", "abx_ddr", "bristol_ddr", "smoke12" ))




intake <- c('ala', 'epa', 'dpa', 'dha','omega3','omega3_noala','omega6')
pfa <- c('epa', 'dpa', 'dha','omega3','omega3_noala','omega6')
maaslin_adjust<-function(adjust){
  dir = paste('./maaslin_results/adjusted_', adjust, '.pcl',sep='')
  Maaslin2(species, 
           totmeta_stn, 
           dir, ##### input
           transform="AST",
           min_abundance=0.0001,
           min_prevalence=0.1,
           random_effects="ID1", 
           fixed_effects = c(paste(adjust, "_pfa", sep = ""),paste(adjust, "_ddr", sep = ""), "agemlvs","calor_fs_dr_ddr", "probx_ddr", "abx_ddr", "bristol_ddr", "smoke12" ))
}

for (a in adjust){
  maaslin_adjust(adjust = a)
}



