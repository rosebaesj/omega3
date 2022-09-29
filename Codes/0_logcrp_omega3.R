# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#+ Proof check on Exposure and Outcome.
#+ Exposure: omega3 intake (without ALA) DDR, FFQ, 10v(cummulative) data are available
#+            need to find wich dietary data to focus.
#+ Outcome: logCRP as the indicator of systematic inflammation
#+ First run brief correlation analysis on 
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

a<- NULL
#### 1. logcrp ~ ddr ####


histogram(meta_bispecies$logcrp_plasma)
histogram(meta_bispecies$omega3_noala_ddr)
histogram(meta_bispecies$logomega3_noala_ddr)

reg <- lmer( logcrp_plasma ~ logomega3_noala_ddr
             +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)

histogram(meta_bispecies$logepa_ddr)
reg <- lmer( logcrp_plasma ~ logepa_ddr
             +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)

histogram(meta_bispecies$logdpa_ddr)
reg <- lmer( logcrp_plasma ~ logdpa_ddr
             +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)

histogram(meta_bispecies$logdha_ddr)
reg <- lmer( logcrp_plasma ~ logdha_ddr
             +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)

histogram(meta_bispecies$logala_ddr)
reg <- lmer( logcrp_plasma ~ logala_ddr
             +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)

histogram(meta_bispecies$logomega3_ddr)
reg <- lmer( logcrp_plasma ~ logomega3_ddr
             +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)



#### 2. logcrp ~ ffq ####

histogram(meta_bispecies$logcrp_plasma)
histogram(meta_bispecies$omega3_noala_ffq)
histogram(meta_bispecies$logomega3_noala_ffq)


meta_bispecies$logomega3_noala_ffq[meta_bispecies$logomega3_noala_ffq%in%c(-Inf,NA, Inf)] <- median(c(meta_bispecies$logomega3_noala_ffq), na.rm = T)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ffq
             +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)

meta_bispecies$logepa_ffq[meta_bispecies$logepa_ffq%in%c(-Inf,NA, Inf)] <- median(c(meta_bispecies$logepa_ffq), na.rm = T)

histogram(meta_bispecies$epa_ffq)
reg <- lmer( logcrp_plasma ~ logepa_ffq
                          +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr  
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)

histogram(meta_bispecies$dpa_ffq)
reg <- lmer( logcrp_plasma ~ logdpa_ffq
                          +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr  
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)

meta_bispecies$logdha_ffq[meta_bispecies$logdha_ffq%in%c(-Inf,NA, Inf)] <- median(c(meta_bispecies$logdha_ffq), na.rm = T)

histogram(meta_bispecies$dha_ffq)
reg <- lmer( logcrp_plasma ~ logdha_ffq
                          +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr  
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)


meta_bispecies$tfat_avg
histogram(meta_bispecies$ala_ffq)
reg <- lmer( logcrp_plasma ~ logala_ffq 
             +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr  
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)

histogram(meta_bispecies$omega3_ffq)
reg <- lmer( logcrp_plasma ~ logomega3_ffq
                          +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr  
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)




#### 3. logcrp ~ 10v ####

histogram(meta_bispecies$logcrp_plasma)
histogram(meta_bispecies$omega3_noala10v)
histogram(meta_bispecies$logomega3_noala10v)

reg <- lmer( logcrp_plasma ~ logomega3_noala10v
             +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)

histogram(meta_bispecies$epa10v)
reg <- lmer( logcrp_plasma ~ logepa10v
             +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr  
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)

histogram(meta_bispecies$dpa10v)
reg <- lmer( logcrp_plasma ~ logdpa10v
             +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr  
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)

meta_bispecies$f



histogram(meta_bispecies$dha10v)
reg <- lmer( logcrp_plasma ~ logdha10v
             +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr  
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)

histogram(meta_bispecies$ala10v)
reg <- lmer( logcrp_plasma ~ logala10v
             +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr  
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)

histogram(meta_bispecies$omega310v)
reg <- lmer( logcrp_plasma ~ logomega310v
             +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr  
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)





#### 4. logcrp ~ pfa ####

histogram(meta_bispecies$logcrp_plasma)
histogram(meta_bispecies$omega3_noala_pfa)
histogram(meta_bispecies$logomega3_noala_pfa)

reg <- lmer( logcrp_plasma ~ logomega3_noala_pfa
             +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)

histogram(meta_bispecies$epa_pfa)
reg <- lmer( logcrp_plasma ~ logepa_pfa
             +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr  
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)

histogram(meta_bispecies$dpa_pfa)
reg <- lmer( logcrp_plasma ~ logdpa_pfa
             +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr  
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)

meta_bispecies$f



histogram(meta_bispecies$dha_pfa)
reg <- lmer( logcrp_plasma ~ logdha_pfa
             +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr  
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)

histogram(meta_bispecies$ala_pfa)
reg <- lmer( logcrp_plasma ~ logala_pfa
             +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr  
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)

histogram(meta_bispecies$omega3_pfa)
reg <- lmer( logcrp_plasma ~ logomega3_pfa
             +fishffq +bmi12+ tfat_avg +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr  
             +(1 | ID1), data=meta_bispecies)
a<- rbind(a, cbind(beta = reg@beta[2],  anova(reg)[1,]))
summary(reg)









a$name <- factor(a$name, levels = rev(a$name), order = T)
a$stars <- cut(a$`Pr(>F)`, breaks=c(-Inf, 0.01, 0.05, 0.10, Inf), label=c("***","**", "*", ""))

write.table(a, file = "./data/logcrp_omega3_ddrffq10v.pcl",  sep = '\t', quote = F, eol = '\n',  row.names=T)


a$data <- factor(c('DDR','DDR','DDR','DDR','DDR','DDR',
               'FFQ','FFQ','FFQ','FFQ','FFQ','FFQ',
               '~2010','~2010','~2010','~2010','~2010','~2010',
               'pfa','pfa','pfa','pfa','pfa','pfa'),
               levels = c('DDR', 'FFQ', '~2010', 'pfa'), order = T)
a$omega3 <- factor(c('logomega3_noala', 'logepa', 'logdpa', 'logdha', 'logala', 'logomega3',
              'logomega3_noala', 'logepa', 'logdpa', 'logdha', 'logala', 'logomega3',
              'logomega3_noala', 'logepa', 'logdpa', 'logdha', 'logala', 'logomega3',
              'logomega3_noala', 'logepa', 'logdpa', 'logdha', 'logala', 'logomega3'),
              levels = rev(c('logomega3_noala', 'logepa', 'logdpa', 'logdha', 'logala', 'logomega3')),
              order = T)



ggplot(a, aes(name, beta)) +
  geom_bar(stat="identity")+
  theme_light()+
  scale_y_continuous()+
  geom_text(aes(label=stars), color="black", size=5) +
  # scale_fill_manual(values = color_code)+ 
  coord_flip()+
  labs(title ="logcrp_logepa_pfa.pcl")+
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

### heatmap #####
ggplot(a, aes(data, omega3)) +
  geom_tile(aes(fill = beta), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=3, show.legend = TRUE) +
  xlab("omega3 intake") +
  labs(title = "Adjusted (bmi, age, calor, abx, probx ...)")+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "right",
        plot.title = element_text(size=10),
        axis.title.x=element_text(size=10,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10, face="italic",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=10,color="black")) +
  labs(fill = "Beta coefficient")


###### scatter plot ######
ggplot(meta_bispecies, aes(logomega3_noala_ddr, logcrp_plasma)) + 
  # ggplot(meta_species, aes(act13m, bmimbs, col=factor(cutx))) + 
  geom_point(size=1, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=1) +
  theme_classic()+
  theme(title = element_blank(),
        plot.title = element_text(face = "bold.italic", size=10),
        legend.position="bottom",
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10))

ggplot(meta_bispecies, aes(logomega3_noala_ffq, logcrp_plasma)) + 
  # ggplot(meta_species, aes(act13m, bmimbs, col=factor(cutx))) + 
  geom_point(size=1, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=1) +
  theme_classic()+
  theme(title = element_blank(),
        plot.title = element_text(face = "bold.italic", size=10),
        legend.position="bottom",
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10))

ggplot(meta_bispecies, aes(logomega3_noala10v, logcrp_plasma)) + 
  # ggplot(meta_species, aes(act13m, bmimbs, col=factor(cutx))) + 
  geom_point(size=1, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=1) +
  theme_classic()+
  theme(title = element_blank(),
        plot.title = element_text(face = "bold.italic", size=10),
        legend.position="bottom",
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10))

ggplot(meta_bispecies, aes(logomega3_noala_pfa, logcrp_plasma)) + 
  # ggplot(meta_species, aes(act13m, bmimbs, col=factor(cutx))) + 
  geom_point(size=1, alpha=0.6, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=1) +
  theme_classic()+
  theme(title = element_blank(),
        plot.title = element_text(face = "bold.italic", size=10),
        legend.position="bottom",
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10))






hist(meta_bispecies$logdpa_pfa)

list <- c('logomega3_noala_ddr', 'logepa_ddr', 'logdpa_ddr', 'logdha_ddr', 'logala_ddr', 'logomega3_ddr',
          'logomega3_noala_ffq','logepa_ffq', 'logdpa_ffq', 'logdha_ffq','logala_ffq','logomega3_ffq',
          'logomega3_noala10v','logepa10v', 'logdpa10v', 'logdha10v','logala10v','logomega310v',
          'logomega3_noala_pfa','logepa_pfa', 'logdpa_pfa', 'logdha_pfa', 'logala_pfa','logomega3_pfa'
)


corr <- NULL
for (l in list){
  spearman<- cor.test(meta_bispecies$logcrp_plasma, meta_bispecies[,l], method = "spearman")
  
  corr <- rbind(corr, c(name = l, rho = spearman$estimate, pval = spearman$p.value))
  
}


corr <- data.frame(corr)
corr$name <- factor(corr$name, levels = rev(corr$name), order = T)
corr$stars <- cut(as.numeric(corr$pval), breaks=c(-Inf, 0.01, 0.05, 0.10, Inf), label=c("***","**", "*", ""))
corr$data <- factor(c('DDR','DDR','DDR','DDR','DDR','DDR',
                   'FFQ','FFQ','FFQ','FFQ','FFQ','FFQ',
                   '~2010','~2010','~2010','~2010','~2010','~2010',
                   'pfa','pfa','pfa','pfa','pfa','pfa'),
                 levels = c('DDR', 'FFQ', '~2010', 'pfa'), order = T)
corr$omega3 <- factor(c('logomega3_noala', 'logepa', 'logdpa', 'logdha', 'logala', 'logomega3',
                     'logomega3_noala', 'logepa', 'logdpa', 'logdha', 'logala', 'logomega3',
                     'logomega3_noala', 'logepa', 'logdpa', 'logdha', 'logala', 'logomega3',
                     'logomega3_noala', 'logepa', 'logdpa', 'logdha', 'logala', 'logomega3'),
                   levels = rev(c('logomega3_noala', 'logepa', 'logdpa', 'logdha', 'logala', 'logomega3')),
                   order = T)



ggplot(corr, aes(data, omega3)) +
  geom_tile(aes(fill = as.numeric(rho.rho)), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=3, show.legend = TRUE) +
  xlab("omega3 intake") +
  labs(title = "Spearman correlation")+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "right",
        plot.title = element_text(size=10),
        axis.title.x=element_text(size=10,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10, face="italic",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=10,color="black")) +
  labs(fill = "Beta coefficient")




#####


list <- c('omega3_noala_ddr', 'epa_ddr', 'dpa_ddr', 'dha_ddr', 'ala_ddr', 'omega3_ddr',
          'omega3_noala_ffq','epa_ffq', 'dpa_ffq', 'dha_ffq','ala_ffq','omega3_ffq',
          'omega3_noala10v','epa10v', 'dpa10v', 'dha10v','ala10v','omega310v',
          'omega3_noala_pfa','epa_pfa', 'dpa_pfa', 'dha_pfa', 'ala_pfa','omega3_pfa'
)


corr <- NULL
for (l in list){
  spearman<- cor.test(meta_bispecies$logcrp_plasma, meta_bispecies[,l], method = "spearman")
  
  corr <- rbind(corr, c(name = l, rho = spearman$estimate, pval = spearman$p.value))
  
}


corr <- data.frame(corr)
corr$name <- factor(corr$name, levels = rev(corr$name), order = T)
corr$stars <- cut(as.numeric(corr$pval), breaks=c(-Inf, 0.01, 0.05, 0.10, Inf), label=c("***","**", "*", ""))
corr$data <- factor(c('DDR','DDR','DDR','DDR','DDR','DDR',
                      'FFQ','FFQ','FFQ','FFQ','FFQ','FFQ',
                      '~2010','~2010','~2010','~2010','~2010','~2010',
                      'pfa','pfa','pfa','pfa','pfa','pfa'),
                    levels = c('DDR', 'FFQ', '~2010', 'pfa'), order = T)
corr$omega3 <- factor(c('omega3_noala', 'epa', 'dpa', 'dha', 'ala', 'omega3',
                        'omega3_noala', 'epa', 'dpa', 'dha', 'ala', 'omega3',
                        'omega3_noala', 'epa', 'dpa', 'dha', 'ala', 'omega3',
                        'omega3_noala', 'epa', 'dpa', 'dha', 'ala', 'omega3'),
                      levels = rev(c('omega3_noala', 'epa', 'dpa', 'dha', 'ala', 'omega3')),
                      order = T)



ggplot(corr, aes(data, omega3)) +
  geom_tile(aes(fill = as.numeric(rho.rho)), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=3, show.legend = TRUE) +
  xlab("omega3 intake") +
  labs(title = "Spearman correlation")+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "right",
        plot.title = element_text(size=10),
        axis.title.x=element_text(size=10,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10, face="italic",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=10,color="black")) +
  labs(fill = "Beta coefficient")





