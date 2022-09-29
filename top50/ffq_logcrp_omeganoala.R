

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
list <- c('logomega3_noala_ffq' , 
          'logepa_pfa', 'logdpa_pfa', 'logdha_pfa', 
          'logAA_pfa', 'logcrp_plasma')
histogram(meta_logspecies$logomega3_noala_ffq)
histogram(meta_logspecies$logcrp_plasma)
histogram(meta_logspecies$AA_pfa)
histogram(meta_logspecies$logepa_pfa)

a <- NULL

#### 1. s__Bacteroides_uniformis ####

histogram(meta_logspecies$s__Bacteroides_uniformis)

table(meta_bispecies$s__Bacteroides_uniformis)

summary(meta_logspecies[meta_bispecies$s__Bacteroides_uniformis==1,]$s__Bacteroides_uniformis)
summary(meta_logspecies[meta_bispecies$s__Bacteroides_uniformis==0,]$s__Bacteroides_uniformis)

reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Bacteroides_uniformis + logomega3_noala_ffq*s__Bacteroides_uniformis 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))
summary(reg)


#### 2. s__Bacteroides_vulgatus ####
histogram(meta_logspecies$s__Bacteroides_vulgatus)
table(meta_bispecies$s__Bacteroides_vulgatus)

summary(meta_logspecies[meta_bispecies$s__Bacteroides_vulgatus==1,]$s__Bacteroides_vulgatus)
summary(meta_logspecies[meta_bispecies$s__Bacteroides_vulgatus==0,]$s__Bacteroides_vulgatus)

reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Bacteroides_vulgatus + logomega3_noala_ffq*s__Bacteroides_vulgatus 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))
summary(reg)



##***3. s__Faecalibacterium_prausnitzii***##

histogram(meta_logspecies$s__Faecalibacterium_prausnitzii)
table(meta_bispecies$s__Faecalibacterium_prausnitzii)

summary(meta_logspecies[meta_bispecies$s__Faecalibacterium_prausnitzii==1,]$s__Faecalibacterium_prausnitzii)
summary(meta_logspecies[meta_bispecies$s__Faecalibacterium_prausnitzii==0,]$s__Faecalibacterium_prausnitzii)

reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Faecalibacterium_prausnitzii + logomega3_noala_ffq*s__Faecalibacterium_prausnitzii 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##4. s__Eubacterium_rectale***####

histogram(meta_logspecies$s__Eubacterium_rectale)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Eubacterium_rectale + logomega3_noala_ffq*s__Eubacterium_rectale 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*5. s__Prevotella_copri***##

histogram(meta_logspecies$s__Prevotella_copri)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Prevotella_copri + logomega3_noala_ffq*s__Prevotella_copri 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*6. s__Alistipes_putredinis***##

histogram(meta_logspecies$s__Alistipes_putredinis)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Alistipes_putredinis + logomega3_noala_ffq*s__Alistipes_putredinis 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*7. s__Bacteroides_dorei***##

histogram(meta_logspecies$s__Bacteroides_dorei)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Bacteroides_dorei + logomega3_noala_ffq*s__Bacteroides_dorei 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


##***8. s__Eubacterium_eligens***##

histogram(meta_logspecies$s__Eubacterium_eligens)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Eubacterium_eligens + logomega3_noala_ffq*s__Eubacterium_eligens 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##9. s__Roseburia_faecis***##

histogram(meta_logspecies$s__Roseburia_faecis)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Roseburia_faecis + logomega3_noala_ffq*s__Roseburia_faecis 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##10. s__Bacteroides_stercoris***##

histogram(meta_logspecies$s__Bacteroides_stercoris)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Bacteroides_stercoris + logomega3_noala_ffq*s__Bacteroides_stercoris 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


##11. s__Eubacterium_siraeum***##

histogram(meta_logspecies$s__Eubacterium_siraeum)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Eubacterium_siraeum + logomega3_noala_ffq*s__Eubacterium_siraeum 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


##12. s__Ruminococcus_bromii***##

histogram(meta_logspecies$s__Ruminococcus_bromii)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Ruminococcus_bromii + logomega3_noala_ffq*s__Ruminococcus_bromii 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


##13. s__Fusicatenibacter_saccharivorans***##

histogram(meta_logspecies$s__Fusicatenibacter_saccharivorans)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Fusicatenibacter_saccharivorans + logomega3_noala_ffq*s__Fusicatenibacter_saccharivorans 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##***14. s__Parabacteroides_distasonis***##

histogram(meta_logspecies$s__Parabacteroides_distasonis)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Parabacteroides_distasonis + logomega3_noala_ffq*s__Parabacteroides_distasonis 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##15. s__Akkermansia_muciniphila***##

histogram(meta_logspecies$s__Akkermansia_muciniphila)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Akkermansia_muciniphila + logomega3_noala_ffq*s__Akkermansia_muciniphila 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##16. s__Bacteroides_eggerthii***##

histogram(meta_logspecies$s__Bacteroides_eggerthii)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Bacteroides_eggerthii + logomega3_noala_ffq*s__Bacteroides_eggerthii 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##17. s__Ruminococcus_bicirculans***##

histogram(meta_logspecies$s__Ruminococcus_bicirculans)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Ruminococcus_bicirculans + logomega3_noala_ffq*s__Ruminococcus_bicirculans 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##18. s__Bacteroides_cellulosilyticus***##

histogram(meta_logspecies$s__Bacteroides_cellulosilyticus)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Bacteroides_cellulosilyticus + logomega3_noala_ffq*s__Bacteroides_cellulosilyticus 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##19. s__Bacteroides_thetaiotaomicron***##

histogram(meta_logspecies$s__Bacteroides_thetaiotaomicron)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Bacteroides_thetaiotaomicron + logomega3_noala_ffq*s__Bacteroides_thetaiotaomicron 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


##20. s__Anaerostipes_hadrus***##

histogram(meta_logspecies$s__Anaerostipes_hadrus)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Anaerostipes_hadrus + logomega3_noala_ffq*s__Anaerostipes_hadrus 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##***21. s__Collinsella_aerofaciens***##

histogram(meta_logspecies$s__Collinsella_aerofaciens)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Collinsella_aerofaciens + logomega3_noala_ffq*s__Collinsella_aerofaciens 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*22. s__Bacteroides_ovatus***##

histogram(meta_logspecies$s__Bacteroides_ovatus)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Bacteroides_ovatus + logomega3_noala_ffq*s__Bacteroides_ovatus 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*23. s__Roseburia_intestinalis***##

histogram(meta_logspecies$s__Roseburia_intestinalis)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Roseburia_intestinalis + logomega3_noala_ffq*s__Roseburia_intestinalis 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##24. s__Alistipes_finegoldii***##

histogram(meta_logspecies$s__Alistipes_finegoldii)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Alistipes_finegoldii + logomega3_noala_ffq*s__Alistipes_finegoldii 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##25. s__Parabacteroides_merdae***##

histogram(meta_logspecies$s__Parabacteroides_merdae)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Parabacteroides_merdae + logomega3_noala_ffq*s__Parabacteroides_merdae 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##***26. s__Barnesiella_intestinihominis***##

histogram(meta_logspecies$s__Barnesiella_intestinihominis)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Barnesiella_intestinihominis + logomega3_noala_ffq*s__Barnesiella_intestinihominis 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##27. s__Bacteroides_caccae***##

histogram(meta_logspecies$s__Bacteroides_caccae)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Bacteroides_caccae + logomega3_noala_ffq*s__Bacteroides_caccae 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##***28. s__Bacteroides_plebeius***##

histogram(meta_logspecies$s__Bacteroides_plebeius)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Bacteroides_plebeius + logomega3_noala_ffq*s__Bacteroides_plebeius 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##*29. s__Dorea_longicatena***##

histogram(meta_logspecies$s__Dorea_longicatena)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Dorea_longicatena + logomega3_noala_ffq*s__Dorea_longicatena 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##30 s__Escherichia_coli***##

histogram(meta_logspecies$s__Escherichia_coli)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Escherichia_coli + logomega3_noala_ffq*s__Escherichia_coli 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



#31. s__Ruminococcus_torques***##

histogram(meta_logspecies$s__Ruminococcus_torques)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Ruminococcus_torques + logomega3_noala_ffq*s__Ruminococcus_torques 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))




##32. s__Lachnospira_pectinoschiza***##

histogram(meta_logspecies$s__Lachnospira_pectinoschiza)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Lachnospira_pectinoschiza + logomega3_noala_ffq*s__Lachnospira_pectinoschiza 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##***33. s__Coprococcus_eutactus***##

histogram(meta_logspecies$s__Coprococcus_eutactus)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Coprococcus_eutactus + logomega3_noala_ffq*s__Coprococcus_eutactus 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##***34. s__Butyrivibrio_crossotus***##

histogram(meta_logspecies$s__Butyrivibrio_crossotus)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Butyrivibrio_crossotus + logomega3_noala_ffq*s__Butyrivibrio_crossotus 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*35. s__Roseburia_inulinivorans***##

histogram(meta_logspecies$s__Roseburia_inulinivorans)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Roseburia_inulinivorans + logomega3_noala_ffq*s__Roseburia_inulinivorans 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*36. s__Bifidobacterium_adolescentis***##

histogram(meta_logspecies$s__Bifidobacterium_adolescentis)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Bifidobacterium_adolescentis + logomega3_noala_ffq*s__Bifidobacterium_adolescentis 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*37. s__Oscillibacter_sp_57_20***##

histogram(meta_logspecies$s__Oscillibacter_sp_57_20)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Oscillibacter_sp_57_20 + logomega3_noala_ffq*s__Oscillibacter_sp_57_20 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##***38. s__Eubacterium_sp_CAG_180***##

histogram(meta_logspecies$s__Eubacterium_sp_CAG_180)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Eubacterium_sp_CAG_180 + logomega3_noala_ffq*s__Eubacterium_sp_CAG_180 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##31. s__Eubacterium_hallii***##

histogram(meta_logspecies$s__Eubacterium_hallii)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Eubacterium_hallii + logomega3_noala_ffq*s__Eubacterium_hallii 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##***39. s__Bacteroides_intestinalis***##

histogram(meta_logspecies$s__Bacteroides_intestinalis)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Bacteroides_intestinalis + logomega3_noala_ffq*s__Bacteroides_intestinalis 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*41. s__Blautia_wexlerae***##

histogram(meta_logspecies$s__Blautia_wexlerae)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Blautia_wexlerae + logomega3_noala_ffq*s__Blautia_wexlerae 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*42. s__Blautia_obeum***##

histogram(meta_logspecies$s__Blautia_obeum)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Blautia_obeum + logomega3_noala_ffq*s__Blautia_obeum 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*43. s__Bifidobacterium_longum***##

histogram(meta_logspecies$s__Bifidobacterium_longum)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Bifidobacterium_longum + logomega3_noala_ffq*s__Bifidobacterium_longum 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##44. s__Coprococcus_comes***##

histogram(meta_logspecies$s__Coprococcus_comes)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Coprococcus_comes + logomega3_noala_ffq*s__Coprococcus_comes 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*45. s__Roseburia_hominis***##

histogram(meta_logspecies$s__Roseburia_hominis)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Roseburia_hominis + logomega3_noala_ffq*s__Roseburia_hominis 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##46. s__Ruminococcus_lactaris***##

histogram(meta_logspecies$s__Ruminococcus_lactaris)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Ruminococcus_lactaris + logomega3_noala_ffq*s__Ruminococcus_lactaris 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##47. s__Bacteroides_fragilis***##

histogram(meta_logspecies$s__Bacteroides_fragilis)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Bacteroides_fragilis + logomega3_noala_ffq*s__Bacteroides_fragilis 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##48. s__Odoribacter_splanchnicus***##

histogram(meta_logspecies$s__Odoribacter_splanchnicus)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Odoribacter_splanchnicus + logomega3_noala_ffq*s__Odoribacter_splanchnicus 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


##*49. s__Butyrivibrio_sp_CAG_318***##

histogram(meta_logspecies$s__Butyrivibrio_sp_CAG_318)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Butyrivibrio_sp_CAG_318 + logomega3_noala_ffq*s__Butyrivibrio_sp_CAG_318 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##50 s__Bacteroides_xylanisolvens***##

histogram(meta_logspecies$s__Bacteroides_xylanisolvens)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq + s__Bacteroides_xylanisolvens + logomega3_noala_ffq*s__Bacteroides_xylanisolvens 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


a$name <- gsub("logomega3_noala_ffq:s__", "", rownames(a))
a$stars <- cut(a$`Pr(>F)`, breaks=c(-Inf, 0.01, 0.05, 0.10, Inf), label=c("***","**", "*", ""))

write.table(a, file = "./data/ffq_int_logcrp_omeganoala.pcl",  sep = '\t', quote = F, eol = '\n',  row.names=T)


ggplot(a, aes(name, beta)) +
  geom_bar(stat="identity")+
  theme_light()+
  scale_y_continuous()+
  geom_text(aes(label=stars), color="black", size=5) +
 # scale_fill_manual(values = color_code)+ 
  coord_flip()+
  labs(title ="logcrp_plasma")+
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








