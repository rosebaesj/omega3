
##############################################################################################################################################
# 1) Purposes: creat heatmap of the assocaitons of measures of physical activity and adiposity with abundances of per microbial species from MaAsLin regression, and prepare data for phylogenetic tree creation in GraPhlan
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
setwd("/udd/nhkwa/mlvs")
options(max.print=1000000)

library(Maaslin2)
library(ggplot2)
library(grid)
library(ggpubr)
library(tidyverse)
library(gridExtra)
library(dplyr)

# read in metadata
# pfa ####
logcrp_noala_pfa <- read.table('./data/interaction_logcrp_o3noala_pfa.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_ala_pfa <- read.table('./data/logcrp_logala_pfa.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_epa_pfa <- read.table('./data/logcrp_logepa_pfa.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_dpa_pfa <- read.table('./data/logcrp_logdpa_pfa.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_dha_pfa <- read.table('./data/logcrp_logdha_pfa.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")


logcrp_noala_pfa$interaction <- "omega3_noala_pfa"
logcrp_ala_pfa$interaction <- "logala_pfa"
logcrp_epa_pfa$interaction <- "logepa_pfa"
logcrp_dpa_pfa$interaction <- "logdpa_pfa"
logcrp_dha_pfa$interaction <- "logdha_pfa"

logcrp_interaction_pfa <- rbind(logcrp_noala_pfa, logcrp_ala_pfa, logcrp_epa_pfa, logcrp_dpa_pfa, logcrp_dha_pfa)


ggplot(logcrp_interaction_pfa, aes(interaction, name)) +
  geom_tile(aes(fill = beta), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=3, show.legend = TRUE) +
  xlab("omega3 intake") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "left",
        plot.title = element_text(size=10),
        axis.title.x=element_text(size=10,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10, face="italic",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=10,color="black")) +
  labs(fill = "Beta coefficient")


#### ddr ####

# read in metadata
logcrp_noala_ddr <- read.table('./data/interaction_logcrp_omeganoala.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_ala_ddr <- read.table('./data/logcrp_logala_ddr.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_epa_ddr <- read.table('./data/logcrp_logepa_ddr.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_dpa_ddr <- read.table('./data/logcrp_logdpa_ddr.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_dha_ddr <- read.table('./data/logcrp_logdha_ddr.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")


logcrp_noala_ddr$interaction <- "omega3_noala_ddr"
logcrp_ala_ddr$interaction <- "logala_ddr"
logcrp_epa_ddr$interaction <- "logepa_ddr"
logcrp_dpa_ddr$interaction <- "logdpa_ddr"
logcrp_dha_ddr$interaction <- "logdha_ddr"

logcrp_interaction_ddr <- rbind(logcrp_noala_ddr, logcrp_ala_ddr, logcrp_epa_ddr, logcrp_dpa_ddr, logcrp_dha_ddr)


ggplot(logcrp_interaction_ddr, aes(interaction, name)) +
  geom_tile(aes(fill = beta), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=3, show.legend = TRUE) +
  xlab("omega3 intake") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "left",
        plot.title = element_text(size=10),
        axis.title.x=element_text(size=10,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10, face="italic",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=10,color="black")) +
  labs(fill = "Beta coefficient")



#### ffq ####

# read in metadata
logcrp_noala_ffq <- read.table('./data/ffq_int_logcrp_omeganoala.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_ala_ffq <- read.table('./data/logcrp_logala_ffq.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_epa_ffq <- read.table('./data/logcrp_logepa_ffq.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_dpa_ffq <- read.table('./data/logcrp_logdpa_ffq.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_dha_ffq <- read.table('./data/logcrp_logdha_ffq.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")


logcrp_noala_ffq$interaction <- "omega3_noala_ffq"
logcrp_ala_ffq$interaction <- "logala_ffq"
logcrp_epa_ffq$interaction <- "logepa_ffq"
logcrp_dpa_ffq$interaction <- "logdpa_ffq"
logcrp_dha_ffq$interaction <- "logdha_ffq"

logcrp_interaction_ffq <- rbind(logcrp_noala_ffq, logcrp_ala_ffq, logcrp_epa_ffq, logcrp_dpa_ffq, logcrp_dha_ffq)


ggplot(logcrp_interaction_ffq, aes(interaction, name)) +
  geom_tile(aes(fill = beta), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=3, show.legend = TRUE) +
  xlab("omega3 intake") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "left",
        plot.title = element_text(size=10),
        axis.title.x=element_text(size=10,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10, face="italic",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=10,color="black")) +
  labs(fill = "Beta coefficient")







#### 10v ####

# read in metadata
logcrp_noala10v <- read.table('./data/logcrp_omeganoala10v.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_ala10v <- read.table('./data/logcrp_logala10v.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_epa10v <- read.table('./data/logcrp_logepa10v.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_dpa10v <- read.table('./data/logcrp_logdpa10v.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_dha10v <- read.table('./data/logcrp_logdha10v.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")


logcrp_noala10v$interaction <- "omega3_noala10v"
logcrp_ala10v$interaction <- "logala10v"
logcrp_epa10v$interaction <- "logepa10v"
logcrp_dpa10v$interaction <- "logdpa10v"
logcrp_dha10v$interaction <- "logdha10v"

logcrp_interaction10v <- rbind(logcrp_noala10v, logcrp_ala10v, logcrp_epa10v, logcrp_dpa10v, logcrp_dha10v)


ggplot(logcrp_interaction10v, aes(interaction, name)) +
  geom_tile(aes(fill = beta), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=3, show.legend = TRUE) +
  xlab("omega3 intake") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "left",
        plot.title = element_text(size=10),
        axis.title.x=element_text(size=10,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10, face="italic",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=10,color="black")) +
  labs(fill = "Beta coefficient")







### total

logcrp_interaction<- rbind(logcrp_interaction_ddr, logcrp_interaction_ffq, logcrp_interaction10v, logcrp_interaction_pfa)

order <- c( "omega3_noala_ddr",
            "omega3_noala_ffq", 
            "omega3_noala10v",
            "omega3_noala_pfa",
            
           "logepa_ddr", 
           "logdpa_ddr",
           "logdha_ddr",
          
           "logepa_ffq",
           "logdpa_ffq",
           "logdha_ffq", 

           "logepa10v",
           "logdpa10v",
           "logdha10v", 
           
           "logepa_pfa",
           "logdpa_pfa",
           "logdha_pfa", 
 
           
           "logala_ddr", 
           "logala_ffq",
           "logala10v",
           "logala_pfa"
           
           )
logcrp_interaction$interaction <- factor(logcrp_interaction$interaction , level = order, order = T )

logcrp_interaction$name <- factor(logcrp_interaction$name, level = rev(gsub("s__","", top50)), order = T )


ggplot(logcrp_interaction, aes(interaction, name)) +
  geom_tile(aes(fill = beta), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=3, show.legend = TRUE) +
  xlab("omega3 intake") +
  labs(title = "int coeff, (to logCRP, pfa adj for intake)")+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "left",
        plot.title = element_text(size=10),
        axis.title.x=element_text(size=10,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10, face="italic",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=10,color="black")) +
  labs(fill = "Beta coefficient")



#### frq ####

# read in metadata
logcrp_noala_frq <- read.table('./data/frq_interaction_logCRP_omeganoala.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_ala_frq <- read.table('./data/frq_int_logCRP_logala_ddr.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_epa_frq <- read.table('./data/frq_int_logCRP_logepa_ddr.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_dpa_frq <- read.table('./data/frq_int_logCRP_logdpa_ddr.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logcrp_dha_frq <- read.table('./data/frq_int_logCRP_logdha_ddr.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")


logcrp_noala_frq$interaction <- "omega3_noala_ddr"
logcrp_ala_frq$interaction <- "logala_ddr"
logcrp_epa_frq$interaction <- "logepa_ddr"
logcrp_dpa_frq$interaction <- "logdpa_ddr"
logcrp_dha_frq$interaction <- "logdha_ddr"

logcrp_interaction_frq <- rbind(logcrp_noala_frq, logcrp_ala_frq, logcrp_epa_frq, logcrp_dpa_frq, logcrp_dha_frq)


ggplot(logcrp_interaction_frq, aes(interaction, name)) +
  geom_tile(aes(fill = beta), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=3, show.legend = TRUE) +
  xlab("omega3 intake") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "left",
        plot.title = element_text(size=10),
        axis.title.x=element_text(size=10,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10, face="italic",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=10,color="black")) +
  labs(fill = "Beta coefficient")



