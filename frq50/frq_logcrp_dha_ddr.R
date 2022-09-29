

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
list <- c('logdha_ddr' , 
          'logepa_pfa', 'logdpa_pfa', 'logdha_pfa', 
          'logcrp_plasma', 'logcrp_plasma')
histogram(meta_logspecies$logdha_ddr)
histogram(meta_logspecies$logcrp_plasma)
histogram(meta_logspecies$AA_pfa)
histogram(meta_logspecies$logepa_pfa)

a <- NULL

#### 1. s__Faecalibacterium_prausnitzii ####
reg <- lmer( logcrp_plasma ~ logdha_ddr + s__Faecalibacterium_prausnitzii + logdha_ddr*s__Faecalibacterium_prausnitzii 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


#### 2. s__Bacteroides_uniformis ####
reg <- lmer( logcrp_plasma ~ logdha_ddr + s__Bacteroides_uniformis + logdha_ddr*s__Bacteroides_uniformis 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##***3. s__Fusicatenibacter_saccharivorans***##
reg <- lmer( logcrp_plasma ~ logdha_ddr + s__Fusicatenibacter_saccharivorans + logdha_ddr*s__Fusicatenibacter_saccharivorans 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


##4. s__Eubacterium_rectale***####

reg <- lmer( logcrp_plasma ~ logdha_ddr + s__Eubacterium_rectale + logdha_ddr*s__Eubacterium_rectale 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*5. s__Anaerostipes_hadrus***##

histogram(meta_logspecies$s__Anaerostipes_hadrus)
reg <- lmer( logcrp_plasma ~ logdha_ddr + s__Anaerostipes_hadrus + logdha_ddr*s__Anaerostipes_hadrus 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##*6. s__Eubacterium_hallii***##

histogram(meta_logspecies$s__Eubacterium_hallii)
reg <- lmer( logcrp_plasma ~ logdha_ddr + s__Eubacterium_hallii + logdha_ddr*s__Eubacterium_hallii 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


##*7. s__Dorea_formicigenerans***##

histogram(meta_logspecies$s__Dorea_formicigenerans)
reg <- lmer( logcrp_plasma ~ logdha_ddr + s__Dorea_formicigenerans + logdha_ddr*s__Dorea_formicigenerans 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


##***8. s__Eubacterium_eligens***##

histogram(meta_logspecies$s__Eubacterium_eligens)
reg <- lmer( logcrp_plasma ~ logdha_ddr + s__Eubacterium_eligens + logdha_ddr*s__Eubacterium_eligens 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##9. s__Agathobaculum_butyriciproducens***##

histogram(meta_logspecies$s__Agathobaculum_butyriciproducens)
reg <- lmer( logcrp_plasma ~ logdha_ddr + s__Agathobaculum_butyriciproducens + logdha_ddr*s__Agathobaculum_butyriciproducens 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##10. s__Roseburia_hominis***##

histogram(meta_logspecies$s__Roseburia_hominis)
reg <- lmer( logcrp_plasma ~ logdha_ddr + s__Roseburia_hominis + logdha_ddr*s__Roseburia_hominis 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


##11. s__Blautia_obeum***##

histogram(meta_logspecies$s__Blautia_obeum)
reg <- lmer( logcrp_plasma ~ logdha_ddr + s__Blautia_obeum + logdha_ddr*s__Blautia_obeum 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


##12. s__Parabacteroides_distasonis***##

histogram(meta_logspecies$s__Parabacteroides_distasonis)
reg <- lmer( logcrp_plasma ~ logdha_ddr + s__Parabacteroides_distasonis + logdha_ddr*s__Parabacteroides_distasonis 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


##13. s__Blautia_wexlerae***##

histogram(meta_logspecies$s__Blautia_wexlerae)
reg <- lmer( logcrp_plasma ~ logdha_ddr + s__Blautia_wexlerae + logdha_ddr*s__Blautia_wexlerae 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##***14. s__Asaccharobacter_celatus***##

histogram(meta_logspecies$s__Asaccharobacter_celatus)
reg <- lmer( logcrp_plasma ~ logdha_ddr + s__Asaccharobacter_celatus + logdha_ddr*s__Asaccharobacter_celatus 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##15. s__Bacteroides_vulgatus***##

histogram(meta_logspecies$s__Bacteroides_vulgatus)
reg <- lmer( logcrp_plasma ~ logdha_ddr + s__Bacteroides_vulgatus + logdha_ddr*s__Bacteroides_vulgatus 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##16. s__Ruthenibacterium_lactatiformans***##

histogram(meta_logspecies$s__Ruthenibacterium_lactatiformans)
reg <- lmer( logcrp_plasma ~ logdha_ddr + s__Ruthenibacterium_lactatiformans + logdha_ddr*s__Ruthenibacterium_lactatiformans 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))



##17. s__Adlercreutzia_equolifaciens***##

histogram(meta_logspecies$s__Adlercreutzia_equolifaciens)
reg <- lmer( logcrp_plasma ~ logdha_ddr + s__Adlercreutzia_equolifaciens + logdha_ddr*s__Adlercreutzia_equolifaciens 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##18. s__Bacteroides_ovatus***##

histogram(meta_logspecies$s__Bacteroides_ovatus)
reg <- lmer( logcrp_plasma ~ logdha_ddr + s__Bacteroides_ovatus + logdha_ddr*s__Bacteroides_ovatus 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






##19. s__Roseburia_inulinivorans***##

histogram(meta_logspecies$s__Roseburia_inulinivorans)
reg <- lmer( logcrp_plasma ~ logdha_ddr + s__Roseburia_inulinivorans + logdha_ddr*s__Roseburia_inulinivorans 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))


##20. s__Dorea_longicatena***##

histogram(meta_logspecies$s__Dorea_longicatena)
reg <- lmer( logcrp_plasma ~ logdha_ddr + s__Dorea_longicatena + logdha_ddr*s__Dorea_longicatena 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),]))






### Draw table ####

a$name <- gsub("logdha_ddr:s__", "", rownames(a))
a$stars <- cut(a$`Pr(>F)`, breaks=c(-Inf, 0.01, 0.05, 0.10, Inf), label=c("***","**", "*", ""))

write.table(a, file = "./data/frq_int_logCRP_logdha_ddr.pcl",  sep = '\t', quote = F, eol = '\n',  row.names=T)


ggplot(a, aes(name, beta)) +
  geom_bar(stat="identity")+
  theme_light()+
  scale_y_continuous()+
  geom_text(aes(label=stars), color="black", size=5) +
  # scale_fill_manual(values = color_code)+ 
  coord_flip()+
  labs(title ="frq_int_logCRP_logdha_ddr.pcl")+
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

