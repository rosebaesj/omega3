

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
library(colourvalues)
library(lme4)
library(lmerTest)


### IMPORT DATA #####
##### + with NAs ####
totmeta_stn <- read.table('./data/totmeta_stn.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
bispecies <- read.table('./data/bispecies_filt.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")

meta_bispecies <- merge(totmeta_stn, bispecies, by = 0)

a <- NULL

#### 1. s__Bacteroides_uniformis ####
reg <- lmer( logepa_pfa ~ logepa_ddr+ s__Bacteroides_uniformis + logepa_ddr*s__Bacteroides_uniformis 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr +calor_fs_dr_ddr
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(beta)],  anova(reg)[nrow(anova(reg)),]))
summary(reg)

hist(meta_bispecies$logepa_pfa)
hist(meta_bispecies$logepa_ddr)
hist(meta_bispecies$epa_ddr)
View(meta_bispecies$logepa_pfa)
hist(meta_bispecies$s__Bacteroides_uniformis)
hist(meta_bispecies$logepa_ddr * meta_bispecies$s__Bacteroides_uniformis)
hist(meta_bispecies$agemlvs)
hist(meta_bispecies$abx_ddr)
hist(meta_bispecies$probx_ddr)
hist(meta_bispecies$bristol_ddr) ## cat or ordinal
hist(meta_bispecies$calor_fs_dr_ddr)



reg <- lmer( logepa_ddr ~  logepa_pfa+ s__Bacteroides_uniformis + logepa_pfa*s__Bacteroides_uniformis 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr +calor_fs_dr_ddr
             +(1 | ID1), data=meta_bispecies)

