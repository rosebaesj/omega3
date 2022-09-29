
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

list <- c('logomega3_noala_ddr' , 
                   'logepa_pfa', 'logdpa_pfa', 'logdha_pfa', 
                  'logAA_pfa', 'logcrp_plasma')

list <- c('logomega3_ddr', 'logala_ddr','logepa_ddr' ,'logdpa_ddr' ,'logdha_ddr', 'logomega3_noala_ddr' , 
  'logala_pfa', 'logepa_pfa', 'logdpa_pfa', 'logdha_pfa', 
  'logAA_pfa', 'logcrp_plasma', 'loghdl_plasma', 'logtc_plasma', 'logtg_plasma')

# scatterplot with x-axis of s__Alistipes_putredinis and y-axis of tdee_pam colored by bmi_dlw
meta_species$omega3_noala_avg

for (i in 1:20){
  for(j in 1:length(list)){
    
    a<- ggplot(meta_species, aes(meta_species[,top20[i]], meta_species[,list[j]])) + 
      geom_point(aes(color = logcrp_plasma),size=2, alpha=0.8)+#, position = position_jitter()) + 
      scale_color_gradientn(colours = c("#fde725ff",  "#dce319ff","#b8de29ff","#29af7fff", "#2d708eff", "#440154ff"))+
      xlab(paste('Relative abundance of ',top20[i], ' (log scale)'))+
      ylab(paste(list[j]))+theme_classic()+
      theme(legend.position="right",
            axis.line = element_line(colour = "black", 
                                     size = 1, linetype = "solid"),
            legend.text = element_text(size=10),
            legend.title = element_blank(),
            axis.text.y=element_text(color = "Black",size=10),
            axis.text.x=element_text(color = "Black",size=10),
            axis.title.y=element_text(color = "Black",size=10),
            axis.title.x=element_text(color = "Black",size=10))
    png(file=paste("figure_generated/scatter_top20/", top20[i], "_", list[j], ".png"),
        width=5,height=4, units="in", res = 1000)
    grid.arrange(a, ncol=1, nrow=1)
    dev.off()
  }
}


for (i in 1:20){
  for(j in 1:length(list)){
    
    a<- 
      ggplot(meta_species, aes(meta_species[,top20[i]], meta_species[,list[j]])) + 
      geom_point(aes(color = logomega3_noala_ddr),size=2, alpha=0.8)+#, position = position_jitter()) + 
      scale_color_gradientn(colours = c("#fde725ff",  "#dce319ff","#b8de29ff","#29af7fff", "#2d708eff", "#440154ff"))+
      xlab(paste('Relative abundance of ',top20[i], ' (log scale)'))+
      ylab(paste(list[j]))+
      theme_classic()+
      labs(title = paste("P_interaction = ", round(p, digits = 5)))+
      theme(plot.title = element_text(color = "Black",size=10, hjust = 0),
            legend.position="right",
            axis.line = element_line(colour = "black", 
                                     size = 1, linetype = "solid"),
            legend.text = element_text(size=10),
            legend.title = element_text(color = "Black",size=10),
            
            axis.text.y=element_text(color = "Black",size=10),
            axis.text.x=element_text(color = "Black",size=10),
            axis.title.y=element_text(color = "Black",size=10),
            axis.title.x=element_text(color = "Black",size=10))
    png(file=paste("figure_generated/scatter_top20_omega3/", top20[i], "_", list[j], ".png"),
        width=5,height=4, units="in", res = 1000)
    grid.arrange(a, ncol=1, nrow=1)
    dev.off()
  }
}

#logcrp ~ bug*omega3_noala-ddr
for (i in 1:20){
  reg<- lmer(formula = logcrp_plasma ~  meta_species[,which( colnames(meta_species)==top20[i] )]*logomega3_noala_ddr +
               agemlvs+calor_fs_dr_ddr+probx_ddr+abx_ddr+bristol_ddr+ smoke12 + (1 | ID1), 
             data=meta_species)
  an<- anova(reg)
  p<- y[nrow(an), ncol(an)]
  
    a<- 
      ggplot(meta_species, aes(meta_species[,top20[i]], meta_species[,list[j]])) + 
      geom_point(aes(color = logcrp_plasma),size=2, alpha=0.8)+#, position = position_jitter()) + 
      scale_color_gradientn(colours = c("#fde725ff",  "#dce319ff","#b8de29ff","#29af7fff", "#2d708eff", "#440154ff"))+
      xlab(paste('Relative abundance of ',top20[i], ' (log scale)'))+
      ylab(paste(list[j]))+
      theme_classic()+
      labs(title = paste("P_interaction = ", p))+
      theme(plot.title = element_text(color = "Black",size=10, hjust = 0),
            legend.position="right",
            axis.line = element_line(colour = "black", 
                                     size = 1, linetype = "solid"),
            legend.text = element_text(size=10),
            legend.title = element_text(color = "Black",size=10),
            
            axis.text.y=element_text(color = "Black",size=10),
            axis.text.x=element_text(color = "Black",size=10),
            axis.title.y=element_text(color = "Black",size=10),
            axis.title.x=element_text(color = "Black",size=10))
    png(file=paste("figure_generated/scatter_top20_P/", top20[i], "_", "logomega3_noala_ddr","_", "logcrp", ".png", sep =""),
        width=5,height=4, units="in", res = 1000)
    grid.arrange(a, ncol=1, nrow=1)
    dev.off()

}


# AA ~ bug*omega

lmer(formula = logAA_pfa ~  meta_bispecies[,top20[i]]+meta_bispecies[,list[j]] +meta_bispecies[,top20[i]]*meta_bispecies[,list[j]]
     +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
     +(1 | ID1), 
     data=meta_bispecies)







meta_species$logAA_pfa
for (j in 1:length(list)){
  j=6
  i=18
  while (i <= 20){
    reg<- lmer(formula = logAA_pfa ~  meta_bispecies[,top20[i]]+meta_bispecies[,list[j]] +meta_bispecies[,top20[i]]*meta_bispecies[,list[j]]
               +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
               +(1 | ID1), 
               data=meta_bispecies)
    an<- anova(reg)
    p<- an[nrow(an), ncol(an)]
    
    a<- 
      ggplot(meta_species, aes(meta_species[,top20[i]], meta_species[,list[j]])) + 
      geom_point(aes(color = logAA_pfa),size=2, alpha=0.8)+#, position = position_jitter()) + 
      scale_color_gradientn(colours = c("#fde725ff",  "#dce319ff","#b8de29ff","#29af7fff", "#2d708eff", "#440154ff"))+
      xlab(paste('Relative abundance of ',top20[i], ' (log scale)'))+
      ylab(paste(list[j]))+
      theme_classic()+
      labs(title = paste("P_interaction = ", p))+
      theme(plot.title = element_text(color = "Black",size=10, hjust = 0),
            legend.position="right",
            axis.line = element_line(colour = "black", 
                                     size = 1, linetype = "solid"),
            legend.text = element_text(size=10),
            legend.title = element_text(color = "Black",size=10),
            
            axis.text.y=element_text(color = "Black",size=10),
            axis.text.x=element_text(color = "Black",size=10),
            axis.title.y=element_text(color = "Black",size=10),
            axis.title.x=element_text(color = "Black",size=10))
    png(file=paste("figure_generated/scatter_top20_logAA_bi/", top20[i], "_", list[j],"_", "logAA_pfa", ".png", sep =""),
        width=5,height=4, units="in", res = 1000)
    grid.arrange(a, ncol=1, nrow=1)
    dev.off()  
    i = i+1
  }
}



meta_species$logAA_pfa
for (j in 1:length(list)){
  j=1
  i=7
  while (i <= 20){
    reg<- lmer(formula = logAA_pfa ~  meta_species[,top20[i]]+meta_species[,list[j]] +meta_species[,top20[i]]*meta_species[,list[j]]
               +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
               +(1 | ID1), 
               data=meta_species)
    an<- anova(reg)
    p<- an[nrow(an), ncol(an)]
    
    a<- 
      ggplot(meta_species, aes(meta_species[,top20[i]], meta_species[,list[j]])) + 
      geom_point(aes(color = logAA_pfa),size=2, alpha=0.8)+#, position = position_jitter()) + 
      scale_color_gradientn(colours = c("#fde725ff",  "#dce319ff","#b8de29ff","#29af7fff", "#2d708eff", "#440154ff"))+
      xlab(paste('Relative abundance of ',top20[i], ' (log scale)'))+
      ylab(paste(list[j]))+
      theme_classic()+
      labs(title = paste("P_interaction = ", p))+
      theme(plot.title = element_text(color = "Black",size=10, hjust = 0),
            legend.position="right",
            axis.line = element_line(colour = "black", 
                                     size = 1, linetype = "solid"),
            legend.text = element_text(size=10),
            legend.title = element_text(color = "Black",size=10),
            
            axis.text.y=element_text(color = "Black",size=10),
            axis.text.x=element_text(color = "Black",size=10),
            axis.title.y=element_text(color = "Black",size=10),
            axis.title.x=element_text(color = "Black",size=10))
    png(file=paste("figure_generated/scatter_top20_logAA/", top20[i], "_", list[j],"_", "logAA_pfa", ".png", sep =""),
        width=5,height=4, units="in", res = 1000)
    grid.arrange(a, ncol=1, nrow=1)
    dev.off()  
    i = i+1
  }
}



ggplot(meta_logspecies, aes(meta_logspecies[,"s__Roseburia_faecis"], meta_logspecies[,"logala_ddr"])) + 
  geom_point(aes(color = logcrp_plasma),size=2, alpha=0.8)+#, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff",  "#dce319ff","#b8de29ff","#29af7fff","#29af7fff", "#2d708eff", "#2d708eff","#440154ff","#440154ff", "#440154ff"))+
  xlab(paste('Relative abundance of s__Roseburia_faecis (log scale)'))+
  ylab(paste("logomega3_noala_pfa"))+
  theme_classic()+
  labs(title = paste("P_interaction = "))+
  theme(plot.title = element_text(color = "Black",size=10, hjust = 0),
        legend.position="right",
        axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"),
        legend.text = element_text(size=10),
        legend.title = element_text(color = "Black",size=10),
        
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10),
        axis.title.y=element_text(color = "Black",size=10),
        axis.title.x=element_text(color = "Black",size=10))





