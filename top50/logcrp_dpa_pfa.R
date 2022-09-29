

rm(list=ls())
getwd()
setwd("./noGit/")
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
species <-read.table(file = './data/species_filt.csv',
                     sep = ',',  row.names = 1,  header = TRUE,   check.names = FALSE)

meta_stn <- read.table('./data/meta_stn.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")

totmeta_stn <- read.table('./data/totmeta_stn.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")




logspecies <-read.table(file = './data/logspecies_filt.pcl', row.names = 1,  header = TRUE, sep='\t',  check.names = FALSE, quote ="")
bispecies <- read.table('./data/bispecies_filt.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")


meta_logspecies <- merge(totmeta_stn, logspecies, by = 0)
meta_bispecies <- merge(totmeta_stn, bispecies, by = 0)

## TOP 20 
order_species <- species [,rev(order(colSums(species)))]
top20 <- colnames(order_species[,1:20])
top50 <- colnames(order_species[,1:50])
order_species[,1:20]

#### WRITE THEM ALL #####
top20
list <- c('logdpa_pfa' , 
          'logepa_pfa', 'logdpa_pfa', 'logdha_pfa', 
          'logAA_pfa', 'logcrp_plasma')
histogram(meta_logspecies$logdpa_pfa)
histogram(meta_logspecies$logcrp_plasma)
histogram(meta_logspecies$AA_pfa)
histogram(meta_logspecies$logepa_pfa)

a <- NULL

#### 1. s__Bacteroides_uniformis ####

histogram(meta_logspecies$s__Bacteroides_uniformis)

table(meta_bispecies$s__Bacteroides_uniformis)

summary(meta_logspecies[meta_bispecies$s__Bacteroides_uniformis==1,]$s__Bacteroides_uniformis)
summary(meta_logspecies[meta_bispecies$s__Bacteroides_uniformis==0,]$s__Bacteroides_uniformis)

reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Bacteroides_uniformis + logdpa_pfa*s__Bacteroides_uniformis 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)

a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))
summary(reg)


#### 2. s__Bacteroides_vulgatus ####
histogram(meta_logspecies$s__Bacteroides_vulgatus)
table(meta_bispecies$s__Bacteroides_vulgatus)

summary(meta_logspecies[meta_bispecies$s__Bacteroides_vulgatus==1,]$s__Bacteroides_vulgatus)
summary(meta_logspecies[meta_bispecies$s__Bacteroides_vulgatus==0,]$s__Bacteroides_vulgatus)

reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Bacteroides_vulgatus + logdpa_pfa*s__Bacteroides_vulgatus 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))
summary(reg)



##***3. s__Faecalibacterium_prausnitzii***##

histogram(meta_logspecies$s__Faecalibacterium_prausnitzii)
table(meta_bispecies$s__Faecalibacterium_prausnitzii)

summary(meta_logspecies[meta_bispecies$s__Faecalibacterium_prausnitzii==1,]$s__Faecalibacterium_prausnitzii)
summary(meta_logspecies[meta_bispecies$s__Faecalibacterium_prausnitzii==0,]$s__Faecalibacterium_prausnitzii)

reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Faecalibacterium_prausnitzii + logdpa_pfa*s__Faecalibacterium_prausnitzii 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##4. s__Eubacterium_rectale***####

histogram(meta_logspecies$s__Eubacterium_rectale)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Eubacterium_rectale + logdpa_pfa*s__Eubacterium_rectale 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*5. s__Prevotella_copri***##

histogram(meta_logspecies$s__Prevotella_copri)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Prevotella_copri + logdpa_pfa*s__Prevotella_copri 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*6. s__Alistipes_putredinis***##

histogram(meta_logspecies$s__Alistipes_putredinis)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Alistipes_putredinis + logdpa_pfa*s__Alistipes_putredinis 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*7. s__Bacteroides_dorei***##

histogram(meta_logspecies$s__Bacteroides_dorei)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Bacteroides_dorei + logdpa_pfa*s__Bacteroides_dorei 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


##***8. s__Eubacterium_eligens***##

histogram(meta_logspecies$s__Eubacterium_eligens)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Eubacterium_eligens + logdpa_pfa*s__Eubacterium_eligens 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##9. s__Roseburia_faecis***##

histogram(meta_logspecies$s__Roseburia_faecis)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Roseburia_faecis + logdpa_pfa*s__Roseburia_faecis 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##10. s__Bacteroides_stercoris***##

histogram(meta_logspecies$s__Bacteroides_stercoris)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Bacteroides_stercoris + logdpa_pfa*s__Bacteroides_stercoris 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


##11. s__Eubacterium_siraeum***##

histogram(meta_logspecies$s__Eubacterium_siraeum)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Eubacterium_siraeum + logdpa_pfa*s__Eubacterium_siraeum 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


##12. s__Ruminococcus_bromii***##

histogram(meta_logspecies$s__Ruminococcus_bromii)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Ruminococcus_bromii + logdpa_pfa*s__Ruminococcus_bromii 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


##13. s__Fusicatenibacter_saccharivorans***##

histogram(meta_logspecies$s__Fusicatenibacter_saccharivorans)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Fusicatenibacter_saccharivorans + logdpa_pfa*s__Fusicatenibacter_saccharivorans 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##***14. s__Parabacteroides_distasonis***##

histogram(meta_logspecies$s__Parabacteroides_distasonis)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Parabacteroides_distasonis + logdpa_pfa*s__Parabacteroides_distasonis 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##15. s__Akkermansia_muciniphila***##

histogram(meta_logspecies$s__Akkermansia_muciniphila)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Akkermansia_muciniphila + logdpa_pfa*s__Akkermansia_muciniphila 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##16. s__Bacteroides_eggerthii***##

histogram(meta_logspecies$s__Bacteroides_eggerthii)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Bacteroides_eggerthii + logdpa_pfa*s__Bacteroides_eggerthii 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##17. s__Ruminococcus_bicirculans***##

histogram(meta_logspecies$s__Ruminococcus_bicirculans)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Ruminococcus_bicirculans + logdpa_pfa*s__Ruminococcus_bicirculans 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##18. s__Bacteroides_cellulosilyticus***##

histogram(meta_logspecies$s__Bacteroides_cellulosilyticus)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Bacteroides_cellulosilyticus + logdpa_pfa*s__Bacteroides_cellulosilyticus 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##19. s__Bacteroides_thetaiotaomicron***##

histogram(meta_logspecies$s__Bacteroides_thetaiotaomicron)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Bacteroides_thetaiotaomicron + logdpa_pfa*s__Bacteroides_thetaiotaomicron 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


##20. s__Anaerostipes_hadrus***##

histogram(meta_logspecies$s__Anaerostipes_hadrus)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Anaerostipes_hadrus + logdpa_pfa*s__Anaerostipes_hadrus 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##***21. s__Collinsella_aerofaciens***##

histogram(meta_logspecies$s__Collinsella_aerofaciens)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Collinsella_aerofaciens + logdpa_pfa*s__Collinsella_aerofaciens 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*22. s__Bacteroides_ovatus***##

histogram(meta_logspecies$s__Bacteroides_ovatus)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Bacteroides_ovatus + logdpa_pfa*s__Bacteroides_ovatus 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*23. s__Roseburia_intestinalis***##

histogram(meta_logspecies$s__Roseburia_intestinalis)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Roseburia_intestinalis + logdpa_pfa*s__Roseburia_intestinalis 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##24. s__Alistipes_finegoldii***##

histogram(meta_logspecies$s__Alistipes_finegoldii)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Alistipes_finegoldii + logdpa_pfa*s__Alistipes_finegoldii 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##25. s__Parabacteroides_merdae***##

histogram(meta_logspecies$s__Parabacteroides_merdae)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Parabacteroides_merdae + logdpa_pfa*s__Parabacteroides_merdae 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##***26. s__Barnesiella_intestinihominis***##

histogram(meta_logspecies$s__Barnesiella_intestinihominis)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Barnesiella_intestinihominis + logdpa_pfa*s__Barnesiella_intestinihominis 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##27. s__Bacteroides_caccae***##

histogram(meta_logspecies$s__Bacteroides_caccae)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Bacteroides_caccae + logdpa_pfa*s__Bacteroides_caccae 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##***28. s__Bacteroides_plebeius***##

histogram(meta_logspecies$s__Bacteroides_plebeius)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Bacteroides_plebeius + logdpa_pfa*s__Bacteroides_plebeius 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##*29. s__Dorea_longicatena***##

histogram(meta_logspecies$s__Dorea_longicatena)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Dorea_longicatena + logdpa_pfa*s__Dorea_longicatena 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##30 s__Escherichia_coli***##

histogram(meta_logspecies$s__Escherichia_coli)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Escherichia_coli + logdpa_pfa*s__Escherichia_coli 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



#31. s__Ruminococcus_torques***##

histogram(meta_logspecies$s__Ruminococcus_torques)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Ruminococcus_torques + logdpa_pfa*s__Ruminococcus_torques 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))




##32. s__Lachnospira_pectinoschiza***##

histogram(meta_logspecies$s__Lachnospira_pectinoschiza)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Lachnospira_pectinoschiza + logdpa_pfa*s__Lachnospira_pectinoschiza 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##***33. s__Coprococcus_eutactus***##

histogram(meta_logspecies$s__Coprococcus_eutactus)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Coprococcus_eutactus + logdpa_pfa*s__Coprococcus_eutactus 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##***34. s__Butyrivibrio_crossotus***##

histogram(meta_logspecies$s__Butyrivibrio_crossotus)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Butyrivibrio_crossotus + logdpa_pfa*s__Butyrivibrio_crossotus 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*35. s__Roseburia_inulinivorans***##

histogram(meta_logspecies$s__Roseburia_inulinivorans)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Roseburia_inulinivorans + logdpa_pfa*s__Roseburia_inulinivorans 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*36. s__Bifidobacterium_adolescentis***##

histogram(meta_logspecies$s__Bifidobacterium_adolescentis)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Bifidobacterium_adolescentis + logdpa_pfa*s__Bifidobacterium_adolescentis 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*37. s__Oscillibacter_sp_57_20***##

histogram(meta_logspecies$s__Oscillibacter_sp_57_20)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Oscillibacter_sp_57_20 + logdpa_pfa*s__Oscillibacter_sp_57_20 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##***38. s__Eubacterium_sp_CAG_180***##

histogram(meta_logspecies$s__Eubacterium_sp_CAG_180)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Eubacterium_sp_CAG_180 + logdpa_pfa*s__Eubacterium_sp_CAG_180 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##31. s__Eubacterium_hallii***##

histogram(meta_logspecies$s__Eubacterium_hallii)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Eubacterium_hallii + logdpa_pfa*s__Eubacterium_hallii 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##***39. s__Bacteroides_intestinalis***##

histogram(meta_logspecies$s__Bacteroides_intestinalis)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Bacteroides_intestinalis + logdpa_pfa*s__Bacteroides_intestinalis 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*41. s__Blautia_wexlerae***##

histogram(meta_logspecies$s__Blautia_wexlerae)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Blautia_wexlerae + logdpa_pfa*s__Blautia_wexlerae 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*42. s__Blautia_obeum***##

histogram(meta_logspecies$s__Blautia_obeum)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Blautia_obeum + logdpa_pfa*s__Blautia_obeum 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*43. s__Bifidobacterium_longum***##

histogram(meta_logspecies$s__Bifidobacterium_longum)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Bifidobacterium_longum + logdpa_pfa*s__Bifidobacterium_longum 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##44. s__Coprococcus_comes***##

histogram(meta_logspecies$s__Coprococcus_comes)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Coprococcus_comes + logdpa_pfa*s__Coprococcus_comes 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*45. s__Roseburia_hominis***##

histogram(meta_logspecies$s__Roseburia_hominis)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Roseburia_hominis + logdpa_pfa*s__Roseburia_hominis 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##46. s__Ruminococcus_lactaris***##

histogram(meta_logspecies$s__Ruminococcus_lactaris)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Ruminococcus_lactaris + logdpa_pfa*s__Ruminococcus_lactaris 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##47. s__Bacteroides_fragilis***##

histogram(meta_logspecies$s__Bacteroides_fragilis)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Bacteroides_fragilis + logdpa_pfa*s__Bacteroides_fragilis 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##48. s__Odoribacter_splanchnicus***##

histogram(meta_logspecies$s__Odoribacter_splanchnicus)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Odoribacter_splanchnicus + logdpa_pfa*s__Odoribacter_splanchnicus 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*49. s__Butyrivibrio_sp_CAG_318***##

histogram(meta_logspecies$s__Butyrivibrio_sp_CAG_318)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Butyrivibrio_sp_CAG_318 + logdpa_pfa*s__Butyrivibrio_sp_CAG_318 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##50 s__Bacteroides_xylanisolvens***##

histogram(meta_logspecies$s__Bacteroides_xylanisolvens)
reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Bacteroides_xylanisolvens + logdpa_pfa*s__Bacteroides_xylanisolvens 
             +logdpa_ddr +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


a$name <- gsub("logdpa_pfa:s__", "", rownames(a))
a$stars <- cut(a$`Pr(>F)`, breaks=c(-Inf, 0.01, 0.05, 0.10, Inf), label=c("***","**", "*", ""))

write.table(a, file = "./data/logcrp_logdpa_pfa.pcl",  sep = '\t', quote = F, eol = '\n',  row.names=T)


ggplot(a, aes(name, beta)) +
  geom_bar(stat="identity")+
  theme_light()+
  scale_y_continuous()+
  geom_text(aes(label=stars), color="black", size=5) +
  # scale_fill_manual(values = color_code)+ 
  coord_flip()+
  labs(title ="logcrp_logdpa_pfa.pcl")+
  theme(axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"),
        legend.position = "none",
        legend.text = element_text(size = 10,color="black"),
        legend.title = element_blank(),
        plot.title = element_text(size=10,color="black"),
        axis.title.x=element_blank(),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10,color="black"),
        axis.text.x = element_text(size=10,color="black"))  








