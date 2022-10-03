
##############################################################################################################################################
# 1) Purposes: conduct interaction analysis between physical activity and abundance of A. putredinis in relation to body weight change
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
setwd("./noGit")
options(max.print=1000000)

library(Maaslin2)
library(ggplot2)
library(grid)
library(ggpubr)
library(tidyverse)
library(gridExtra)
library(colourvalues)

# read in meta data
totmeta_stn <- read.table('./data/totmeta_stn.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")

# read in taxonomy data
species <-read.table(file = './data/species_filt.csv',
                     sep = ',',  row.names = 1,  header = TRUE,   check.names = FALSE)

meta_species <- merge(totmeta_stn, species, by = 0)
logspecies <-read.table(file = './data/logspecies_filt.pcl', row.names = 1,  header = TRUE, sep='\t',  check.names = FALSE, quote ="")
bispecies <- read.table('./data/bispecies_filt.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")


meta_logspecies <- merge(totmeta_stn, logspecies, by = 0)
meta_bispecies <- merge(totmeta_stn, bispecies, by = 0)


top20<- c(
"s__Bacteroides_uniformis"           ,"s__Bacteroides_vulgatus"     ,       "s__Faecalibacterium_prausnitzii",    "s__Eubacterium_rectale"  ,          
 "s__Prevotella_copri"              ,"s__Alistipes_putredinis"       ,     "s__Bacteroides_dorei"             ,  "s__Eubacterium_eligens"    ,        
 "s__Roseburia_faecis"              ,  "s__Bacteroides_stercoris"     ,      "s__Eubacterium_siraeum"          ,   "s__Ruminococcus_bromii"   ,         
 "s__Fusicatenibacter_saccharivorans", "s__Parabacteroides_distasonis" ,     "s__Akkermansia_muciniphila"       ,  "s__Bacteroides_eggerthii"  ,        
 "s__Ruminococcus_bicirculans"        ,"s__Bacteroides_cellulosilyticus",    "s__Bacteroides_thetaiotaomicron"   , "s__Anaerostipes_hadrus" 
)

frq20 <- c(
  "s__Faecalibacterium_prausnitzii"  ,  "s__Bacteroides_uniformis",           "s__Fusicatenibacter_saccharivorans" ,"s__Eubacterium_rectale"     ,       
   "s__Anaerostipes_hadrus"           ,  "s__Eubacterium_hallii"   ,           "s__Dorea_formicigenerans" ,          "s__Eubacterium_eligens"     ,       
  "s__Agathobaculum_butyriciproducens" ,"s__Roseburia_hominis"      ,         "s__Blautia_obeum"           ,        "s__Parabacteroides_distasonis",     
 "s__Blautia_wexlerae"            ,    "s__Asaccharobacter_celatus"  ,       "s__Bacteroides_vulgatus"      ,      "s__Ruthenibacterium_lactatiformans",
   "s__Adlercreutzia_equolifaciens",     "s__Bacteroides_ovatus"      ,        "s__Roseburia_inulinivorans"  ,       "s__Dorea_longicatena" 
)






# scatterplot with x-axis of s__Alistipes_putredinis and y-axis of tdee_pam colored by bmi_dlw
i=1

i<- i+1
bug <- which(colnames(meta_logspecies)==top20[i])

ggplot(meta_logspecies, aes(meta_logspecies[,which(colnames(meta_logspecies)==top20[i])], logomega3_noala_ddr)) + 
  geom_point(aes(color = logcrp_plasma),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff", "#b8de29ff","#29af7fff", "#2d708eff", "#440154ff", "#440154ff"))+
  xlab(top20[i])+
  #ylab("Recent total PA (MET-hours/week)")+
  theme_classic()+
  theme(legend.position="none",
        axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10),
        axis.title.y=element_text(color = "Black",size=10),
        axis.title.x=element_text(color = "Black",size=10))

ggplot(meta_species, aes(meta_species[,which(colnames(meta_species)==top20[i])], logomega3_noala_ddr)) + 
  geom_point(aes(color = logcrp_plasma),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff",  "#b8de29ff","#29af7fff", "#2d708eff", "#440154ff", "#440154ff"))+
  xlab(top20[i])+
  #ylab("Recent total PA (MET-hours/week)")+
  theme_classic()+
  theme(legend.position="none",
        axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10),
        axis.title.y=element_text(color = "Black",size=10),
        axis.title.x=element_text(color = "Black",size=10))

i=1

i<- i+1
bug <- which(colnames(meta_logspecies)==frq20[i])

ggplot(meta_logspecies, aes(meta_logspecies[,which(colnames(meta_logspecies)==frq20[i])], logomega3_noala_ddr)) + 
  geom_point(aes(color = logcrp_plasma),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff", "#b8de29ff","#29af7fff", "#2d708eff", "#440154ff", "#440154ff"))+
  xlab(frq20[i])+
  #ylab("Recent total PA (MET-hours/week)")+
  theme_classic()+
  theme(legend.position="none",
        axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10),
        axis.title.y=element_text(color = "Black",size=10),
        axis.title.x=element_text(color = "Black",size=10))

ggplot(meta_species, aes(meta_species[,which(colnames(meta_species)==frq20[i])], logomega3_noala_ddr)) + 
  geom_point(aes(color = logcrp_plasma),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff",  "#b8de29ff","#29af7fff", "#2d708eff", "#440154ff", "#440154ff"))+
  xlab(frq20[i])+
  #ylab("Recent total PA (MET-hours/week)")+
  theme_classic()+
  theme(legend.position="none",
        axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10),
        axis.title.y=element_text(color = "Black",size=10),
        axis.title.x=element_text(color = "Black",size=10))


#############
logcrp_o3noala_ddr <- read.table('./data/interaction_logcrp_omeganoala.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_o3noala_ddr[i,'Pr..F.']
i<- 0 

{
i<- i+1
bug <- which(colnames(meta_logspecies)==top20[i])
ggplot(meta_bispecies, aes(logomega3_noala_ddr, logcrp_plasma, col=factor(meta_bispecies[,which(colnames(meta_bispecies)==top20[i])]))) + 
  geom_point(size=2, alpha=0.4, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=1) +
  ggtitle(top20[i])+
  labs(subtitle = paste( logcrp_o3noala_ddr[i,'stars'], "P_int = ", round (logcrp_o3noala_ddr[i,'Pr..F.'], 5)))+
  #xlab("Physical activity level (MET-hours/week)")+
  #ylab(expression(BMI~at~stool~collection~(kg/(m^2))))+
  theme_classic()+
  theme(title = element_blank(),
        plot.title = element_text(face = "bold.italic", size=10),
        plot.subtitle = element_text(size=10),
        legend.position="none",
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10))
}


i<- 0 

{
  i<- i+1
  bug <- which(colnames(meta_logspecies)==frq20[i])
  ggplot(meta_bispecies, aes(logomega3_noala_ddr, logcrp_plasma, col=factor(meta_bispecies[,which(colnames(meta_bispecies)==frq20[i])]))) + 
    geom_point(size=2, alpha=0.4, position = position_jitter()) +
    geom_smooth(method="lm", se=TRUE, size=1) +
    ggtitle(frq20[i])+
    labs(subtitle = paste( logcrp_o3noala_ddr[i,'stars'], "P_int = ", round (logcrp_o3noala_ddr[i,'Pr..F.'], 5)))+
    #xlab("Physical activity level (MET-hours/week)")+
    #ylab(expression(BMI~at~stool~collection~(kg/(m^2))))+
    theme_classic()+
    theme(title = element_blank(),
          plot.title = element_text(face = "bold.italic", size=10),
          plot.subtitle = element_text(size=10),
          legend.position="none",
          legend.text = element_text(size=10),
          legend.title = element_blank(),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10))
}


####

#############
reg <- lmer( logcrp_plasma ~ logomega3_noala_ddr + s__Faecalibacterium_prausnitzii + logomega3_noala_ddr*s__Faecalibacterium_prausnitzii 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr +bmi12
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
anova(reg)[nrow(anova(reg)),ncol(anova(reg))]
meta_bispecies$s__
ggplot(meta_bispecies, aes(logomega3_noala_ddr, logAA_pfa, col=factor(s__Faecalibacterium_prausnitzii))) + 
  geom_point(size=2, alpha=0.4, position = position_jitter()) +
  geom_smooth(method="lm", se=TRUE, size=1) +
  ggtitle("s__Faecalibacterium_prausnitzii")+
  labs( subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
  theme_classic()+
  theme(title = element_blank(),
        plot.title = element_text(face = "bold.italic", size=10),
        plot.subtitle = element_text(size=10),
        legend.position="none",
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10))




reg <- lmer( logAA_pfa ~ logomega3_noala_ddr + s__Faecalibacterium_prausnitzii + logomega3_noala_ddr*s__Faecalibacterium_prausnitzii 
             #+agemlvs#+abx_ddr#+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_bispecies)
summary(reg)
a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
anova(reg)[nrow(anova(reg)),ncol(anova(reg))]
meta_bispecies$s__
  ggplot(meta_bispecies, aes(logomega3_noala_ddr, logAA_pfa, col=factor(s__Faecalibacterium_prausnitzii))) + 
    geom_point(size=2, alpha=0.4, position = position_jitter()) +
    geom_smooth(method="lm", se=TRUE, size=1) +
    ggtitle("s__Faecalibacterium_prausnitzii")+
    labs( subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
    theme_classic()+
    theme(title = element_blank(),
          plot.title = element_text(face = "bold.italic", size=10),
          plot.subtitle = element_text(size=10),
          legend.position="none",
          legend.text = element_text(size=10),
          legend.title = element_blank(),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10))

  reg <- lmer( AA_pfa ~ logomega3_noala_ddr + s__Faecalibacterium_prausnitzii + logomega3_noala_ddr*s__Faecalibacterium_prausnitzii 
               #+agemlvs#+abx_ddr#+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
               +(1 | ID1), data=meta_bispecies)
  summary(reg)
  a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
  anova(reg)[nrow(anova(reg)),ncol(anova(reg))]
  meta_bispecies$s__
  ggplot(meta_bispecies, aes(logomega3_noala_ddr, AA_pfa, col=factor(s__Faecalibacterium_prausnitzii))) + 
    geom_point(size=2, alpha=0.4, position = position_jitter()) +
    geom_smooth(method="lm", se=TRUE, size=1) +
    ggtitle("s__Faecalibacterium_prausnitzii")+
    labs( subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
    theme_classic()+
    theme(title = element_blank(),
          plot.title = element_text(face = "bold.italic", size=10),
          plot.subtitle = element_text(size=10),
          legend.position="none",
          legend.text = element_text(size=10),
          legend.title = element_blank(),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10))
  
  reg <- lmer( logomega3_noala_pfa ~ logomega3_noala_ddr + s__Faecalibacterium_prausnitzii + logomega3_noala_ddr*s__Faecalibacterium_prausnitzii 
               #+agemlvs#+abx_ddr#+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
               +(1 | ID1), data=meta_bispecies)
  summary(reg)
  a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
  anova(reg)[nrow(anova(reg)),ncol(anova(reg))]
  meta_bispecies$s__
  ggplot(meta_bispecies, aes(logomega3_noala_ddr, logomega3_noala_pfa, col=factor(s__Faecalibacterium_prausnitzii))) + 
    geom_point(size=2, alpha=0.4, position = position_jitter()) +
    geom_smooth(method="lm", se=TRUE, size=1) +
    ggtitle("s__Faecalibacterium_prausnitzii")+
    labs( subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
    theme_classic()+
    theme(title = element_blank(),
          plot.title = element_text(face = "bold.italic", size=10),
          plot.subtitle = element_text(size=10),
          legend.position="none",
          legend.text = element_text(size=10),
          legend.title = element_blank(),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10))
  
  
  
  #AA/EPA
  
  meta_bispecies$AAtoEPA <- meta_bispecies$AA_pfa / meta_bispecies$epa_pfa
  meta_bispecies$logAAtoEPA <- log(meta_bispecies$AAtoEPA)
  
  
  reg <- lmer( logAAtoEPA ~ logomega3_noala_ddr + s__Faecalibacterium_prausnitzii + logomega3_noala_ddr*s__Faecalibacterium_prausnitzii 
               #+agemlvs#+abx_ddr#+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
               +(1 | ID1), data=meta_bispecies)
  summary(reg)
  a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
  anova(reg)[nrow(anova(reg)),ncol(anova(reg))]
  meta_bispecies$s__
  ggplot(meta_bispecies, aes(logomega3_noala_ddr, logAAtoEPA, col=factor(s__Faecalibacterium_prausnitzii))) + 
    geom_point(size=2, alpha=0.4, position = position_jitter()) +
    geom_smooth(method="lm", se=TRUE, size=1) +
    ggtitle("s__Faecalibacterium_prausnitzii")+
    labs( subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
    theme_classic()+
    theme(title = element_blank(),
          plot.title = element_text(face = "bold.italic", size=10),
          plot.subtitle = element_text(size=10),
          legend.position="none",
          legend.text = element_text(size=10),
          legend.title = element_blank(),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10))
  
  
  
  
  #HDL
  
  
  reg <- lmer( loghdl_plasma ~ logomega3_noala_ddr + s__Faecalibacterium_prausnitzii + logomega3_noala_ddr*s__Faecalibacterium_prausnitzii 
               +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
               +(1 | ID1), data=meta_bispecies)
  summary(reg)
  a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
  anova(reg)[nrow(anova(reg)),ncol(anova(reg))]
  meta_bispecies$s__
  ggplot(meta_bispecies, aes(logomega3_noala_ddr, loghdl_plasma, col=factor(s__Faecalibacterium_prausnitzii))) + 
    geom_point(size=2, alpha=0.4, position = position_jitter()) +
    geom_smooth(method="lm", se=TRUE, size=1) +
    ggtitle("s__Faecalibacterium_prausnitzii")+
    labs( subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
    theme_classic()+
    theme(title = element_blank(),
          plot.title = element_text(face = "bold.italic", size=10),
          plot.subtitle = element_text(size=10),
          legend.position="none",
          legend.text = element_text(size=10),
          legend.title = element_blank(),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10))
  
  
  
  # tc
  reg <- lmer( logtg_plasma ~ logomega3_noala_ddr + s__Faecalibacterium_prausnitzii + logomega3_noala_ddr*s__Faecalibacterium_prausnitzii 
               +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
               +(1 | ID1), data=meta_bispecies)
  summary(reg)
  a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
  anova(reg)[nrow(anova(reg)),ncol(anova(reg))]
  meta_bispecies$s__
  ggplot(meta_bispecies, aes(logomega3_noala_ddr, logtg_plasma, col=factor(s__Faecalibacterium_prausnitzii))) + 
    geom_point(size=2, alpha=0.4, position = position_jitter()) +
    geom_smooth(method="lm", se=TRUE, size=1) +
    ggtitle("s__Faecalibacterium_prausnitzii")+
    labs( subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
    theme_classic()+
    theme(title = element_blank(),
          plot.title = element_text(face = "bold.italic", size=10),
          plot.subtitle = element_text(size=10),
          legend.position="none",
          legend.text = element_text(size=10),
          legend.title = element_blank(),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10))

  
  
  # dpa
  reg <- lmer( logala_pfa ~ logomega3_noala_ddr + s__Faecalibacterium_prausnitzii + logomega3_noala_ddr*s__Faecalibacterium_prausnitzii 
               #+agemlvs#+abx_ddr+#probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
               +(1 | ID1), data=meta_bispecies)
  summary(reg)
  a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
  anova(reg)[nrow(anova(reg)),ncol(anova(reg))]
  meta_bispecies$s__
  ggplot(meta_bispecies, aes(logomega3_noala_ddr, logala_pfa, col=factor(s__Faecalibacterium_prausnitzii))) + 
    geom_point(size=2, alpha=0.4, position = position_jitter()) +
    geom_smooth(method="lm", se=TRUE, size=1) +
    ggtitle("s__Faecalibacterium_prausnitzii")+
    labs( subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
    theme_classic()+
    theme(title = element_blank(),
          plot.title = element_text(face = "bold.italic", size=10),
          plot.subtitle = element_text(size=10),
          legend.position="none",
          legend.text = element_text(size=10),
          legend.title = element_blank(),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10))

  
  
  # CRP~ pfa
  reg <- lmer( logcrp_plasma ~ logala_pfa + s__Faecalibacterium_prausnitzii + logala_pfa*s__Faecalibacterium_prausnitzii 
               +bristol_ddr+ agemlvs+abx_ddr+probx_ddr+smoke12 +calor_fs_dr_ddr 
               +(1 | ID1), data=meta_bispecies)
  summary(reg)
  a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
  anova(reg)[nrow(anova(reg)),ncol(anova(reg))]
  meta_bispecies$s__
  ggplot(meta_bispecies, aes(logala_pfa, logcrp_plasma, col=factor(s__Faecalibacterium_prausnitzii))) + 
    geom_point(size=2, alpha=0.4, position = position_jitter()) +
    geom_smooth(method="lm", se=TRUE, size=1) +
    ggtitle("s__Faecalibacterium_prausnitzii")+
    labs( subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
    theme_classic()+
    theme(title = element_blank(),
          plot.title = element_text(face = "bold.italic", size=10),
          plot.subtitle = element_text(size=10),
          legend.position="none",
          legend.text = element_text(size=10),
          legend.title = element_blank(),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10))
  
  reg <- lmer( logcrp_plasma ~ logepa_pfa + s__Faecalibacterium_prausnitzii + logepa_pfa*s__Faecalibacterium_prausnitzii 
               +bristol_ddr+ agemlvs+abx_ddr+probx_ddr+smoke12 +calor_fs_dr_ddr 
               +(1 | ID1), data=meta_bispecies)
  summary(reg)
  a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
  anova(reg)[nrow(anova(reg)),ncol(anova(reg))]
  meta_bispecies$s__
  ggplot(meta_bispecies, aes(logepa_pfa, logcrp_plasma, col=factor(s__Faecalibacterium_prausnitzii))) + 
    geom_point(size=2, alpha=0.4, position = position_jitter()) +
    geom_smooth(method="lm", se=TRUE, size=1) +
    ggtitle("s__Faecalibacterium_prausnitzii")+
    labs( subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
    theme_classic()+
    theme(title = element_blank(),
          plot.title = element_text(face = "bold.italic", size=10),
          plot.subtitle = element_text(size=10),
          legend.position="none",
          legend.text = element_text(size=10),
          legend.title = element_blank(),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10))
  
  reg <- lmer( logcrp_plasma ~ logdpa_pfa + s__Faecalibacterium_prausnitzii + logdpa_pfa*s__Faecalibacterium_prausnitzii 
               +bristol_ddr+ agemlvs+abx_ddr+probx_ddr+smoke12 +calor_fs_dr_ddr 
               +(1 | ID1), data=meta_bispecies)
  summary(reg)
  a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
  anova(reg)[nrow(anova(reg)),ncol(anova(reg))]
  meta_bispecies$s__
  ggplot(meta_bispecies, aes(logdpa_pfa, logcrp_plasma, col=factor(s__Faecalibacterium_prausnitzii))) + 
    geom_point(size=2, alpha=0.4, position = position_jitter()) +
    geom_smooth(method="lm", se=TRUE, size=1) +
    ggtitle("s__Faecalibacterium_prausnitzii")+
    labs( subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
    theme_classic()+
    theme(title = element_blank(),
          plot.title = element_text(face = "bold.italic", size=10),
          plot.subtitle = element_text(size=10),
          legend.position="none",
          legend.text = element_text(size=10),
          legend.title = element_blank(),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10))

  reg <- lmer( logcrp_plasma ~ logdha_pfa + s__Faecalibacterium_prausnitzii + logdha_pfa*s__Faecalibacterium_prausnitzii 
               +bristol_ddr+ agemlvs+abx_ddr+probx_ddr+smoke12 +calor_fs_dr_ddr 
               +(1 | ID1), data=meta_bispecies)
  summary(reg)
  a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
  anova(reg)[nrow(anova(reg)),ncol(anova(reg))]
  meta_bispecies$s__
  ggplot(meta_bispecies, aes(logdha_pfa, logcrp_plasma, col=factor(s__Faecalibacterium_prausnitzii))) + 
    geom_point(size=2, alpha=0.4, position = position_jitter()) +
    geom_smooth(method="lm", se=TRUE, size=1) +
    ggtitle("s__Faecalibacterium_prausnitzii")+
    labs( subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
    theme_classic()+
    theme(title = element_blank(),
          plot.title = element_text(face = "bold.italic", size=10),
          plot.subtitle = element_text(size=10),
          legend.position="none",
          legend.text = element_text(size=10),
          legend.title = element_blank(),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10))
  
  meta_bispecies$bmi12
  reg <- lmer( logdha_pfa ~ logdha_ddr+ (1 | ID1), data = meta_bispecies)
  
  reg <- lmer( bmi12 ~ logomega3_pfa #+ s__Faecalibacterium_prausnitzii + logdha_pfa*s__Faecalibacterium_prausnitzii 
               #+bristol_ddr+ agemlvs+abx_ddr+probx_ddr+smoke12 +calor_fs_dr_ddr +
               , data=meta_bispecies)
  summary(reg)
  a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
  anova(reg)
  meta_bispecies$s__
  ggplot(meta_bispecies, aes(logomega3_pfa, bmi12, col=factor(s__Faecalibacterium_prausnitzii))) + 
    geom_point(size=2, alpha=0.4, position = position_jitter()) +
    geom_smooth(method="lm", se=TRUE, size=1) +
    ggtitle("s__Faecalibacterium_prausnitzii")+
    labs( subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
    theme_classic()+
    theme(title = element_blank(),
          plot.title = element_text(face = "bold.italic", size=10),
          plot.subtitle = element_text(size=10),
          legend.position="none",
          legend.text = element_text(size=10),
          legend.title = element_blank(),
          axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10))
  
  