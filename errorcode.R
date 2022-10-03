

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
reg <- lmer( epa_pfa ~ logepa_ddr+ s__Bacteroides_uniformis + logepa_ddr*s__Bacteroides_uniformis 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr +calor_fs_dr_ddr
             +(1 | ID1), data=meta_bispecies)

reg <- glmer( logepa_pfa ~ logepa_ddr+ s__Bacteroides_uniformis + logepa_ddr*s__Bacteroides_uniformis 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr +calor_fs_dr_ddr
             +(1 | ID1), data=meta_bispecies)


od_glmer = glmer(y1 ~ procedure + (1|id), data = data_g2, family = binomial)
summary(mod_glmer)

glmb = bglmer(y1 ~ procedure + (1|id), data=data_g2, family=binomial,
              fixef.prior = normal(cov = diag(9,2)))

pairs(emmeans(glmb, ~ procedure))

a<- rbind(a, cbind(beta = reg@beta[length(beta)],  anova(reg)[nrow(anova(reg)),]))
summary(reg)
meta_bispecies$epa_pfa
Bacter0 <- meta_bispecies %>%
  filter(s__Bacteroides_uniformis ==0)

Bacter1 <- meta_bispecies %>%
  filter(s__Bacteroides_uniformis ==1)

hist(meta_bispecies$epa_pfa)
hist(Bacter1$logepa_pfa)

hist(meta_bispecies$logepa_pfa)
hist(meta_bispecies$logepa_ddr)
hist(meta_bispecies$epa_ddr)
View(meta_bispecies$logepa_pfa)
hist(meta_bispecies$s__Bacteroides_uniformis)
hist(meta_bispecies$epa_ddr * meta_bispecies$s__Bacteroides_uniformis)
hist(meta_bispecies$agemlvs)
hist(meta_bispecies$abx_ddr)
hist(meta_bispecies$probx_ddr)
hist(meta_bispecies$bristol_ddr) ## cat or ordinal
hist(meta_bispecies$calor_fs_dr_ddr)



reg <- lmer( logepa_ddr ~  logepa_pfa+ s__Bacteroides_uniformis + logepa_pfa*s__Bacteroides_uniformis 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr +calor_fs_dr_ddr
             +(1 | ID1), data=meta_bispecies)

