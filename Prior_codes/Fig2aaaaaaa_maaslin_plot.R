
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
# read in metadata
totmeta_stn <- read.table('./data/totmeta_stn.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")

# read in taxonomy data
species <-read.table(file = './data/species_filt.csv',
                     sep = ',',  row.names = 1,  header = TRUE,   check.names = FALSE)

# combine results form maasline regressions across different dietary variables

exposures <- c('omega3_ddr', 'ala_ddr','epa_ddr' ,'dpa_ddr' ,'dha_ddr', 'omega3_noala_ddr' , 
              'omega310v', 'ala10v', 'epa10v', 'dpa10v', 'dha10v', 'omega3_noala10v')

exposures_ffq <- c('omega3_ddr', 'ala_ddr','epa_ddr' ,'dpa_ddr' ,'dha_ddr', 'omega3_noala_ddr' , 
                   'omega3_ffq', 'ala_ffq','epa_ffq' ,'dpa_ffq' ,'dha_ffq', 'omega3_noala_ffq' , 
               'omega310v', 'ala10v', 'epa10v', 'dpa10v', 'dha10v', 'omega3_noala10v')




logexposures <- paste('log', exposures, sep = "")

outcomes <- c('logcrp_plasma', 'loghdl_plasma', 'logtc_plasma', 'logtg_plasma', 
              'AA_pfa', 'ala_pfa', 'epa_pfa', 'dpa_pfa', 'dha_pfa', 'omega3_pfa',
              'logAA_pfa', 'logala_pfa', 'logepa_pfa', 'logdpa_pfa', 'logdha_pfa', 'logomega3_noala_pfa')

adjust <- c('ala', 'epa', 'dpa', 'dha','omega3','omega3_noala','omega6')

### get significant results ###

sig_results<- function(exposure){
  dir = paste('./maaslin_results/', exposure, '.pcl/significant_results.tsv',sep='')
  result<-read.table(file = dir, sep = '\t', header = TRUE, check.names=FALSE)
  sig_result <- subset(result, metadata==exposure, select =c(feature, metadata))
  return(sig_result)
}


#### all results ####

all_results <- function(exposure, label){
  dir = paste('./maaslin_results/', exposure, '.pcl/all_results.tsv',sep='')
  result<-read.table(file =dir, sep = '\t', header = TRUE, check.names=FALSE)
  all_result <- subset(result, metadata==exposure, select =-c(value, N))
  all_result$meta <- label
  return(all_result)
}


# read in table of species matched to phyla
all_species <- read_tsv("data/metaphlan_taxonomic_profiles.tsv")
names(all_species)[names(all_species) == '# taxonomy'] <- 'Sample'

all_species <-all_species %>%
  separate(Sample, c("kingdom","phylum","class" ,"order","family","genus" ,"species" ), 
           sep = '\\|', remove = TRUE)

all_species_name <- subset(all_species,!is.na(species))[1:7]
all_species_name <- as.data.frame(all_species_name) #weird formatting change to maintain rows
rownames(all_species_name)<-all_species_name$species

colnames(all_species_name)[which(colnames(all_species_name) == 'species')] <- 'feature'

# vector <- all_species_name %>%
#   arrange(feature) %>% arrange(genus) %>% arrange(family) %>% arrange(order) %>% arrange(class) %>% arrange(phylum) %>% arrange(kingdom)
#   
# vector$feature <- substring(all_species_name$feature, 4)
# vector$phylum <- substring(all_species_name$phylum, 4)



####### feature by feature #####
plot_maaslin <- function(listofmeta, by = by){
# listofmeta <- outcomes
# by <- 'logcrp_plasma'
  join <- sig_results(listofmeta[1])
  for (i in 2:length(listofmeta)) {
    join<-full_join(join, sig_results(listofmeta[i]), by="feature")
  }
  sigf_exp<-subset(join, select=feature)
  
  
  bind <- left_join(sigf_exp, all_results(listofmeta[1], label = paste(listofmeta[1], "l", sep = "")), by="feature")
  for (i in 2:length(listofmeta)) {
    bind <- rbind(bind, left_join(sigf_exp, all_results(listofmeta[i], label = paste(listofmeta[i], "l", sep = "")), by="feature"))
  }
  
  bind<-left_join(bind, all_species_name, by="feature")
  
  bind$feature <- substring(bind$feature, 4)
  bind$phylum <- substring(bind$phylum, 4)
  
  bind<-bind[order(bind$meta, bind$phylum, bind$feature),]
  vector <- bind %>%
    arrange(feature) %>% arrange(genus) %>% arrange(family) %>% arrange(order) %>% arrange(class) %>% arrange(phylum) %>% arrange(kingdom) %>% arrange(meta) %>%
    slice(1:(nrow(bind)/length(listofmeta)))
  
  rownames(bind) <- 1:nrow(bind)

  bind$stars <- cut(bind$qval, breaks=c(-Inf, 0.01, 0.05, 0.1, 0.25, Inf), label=c("****", "***", "**", "*", ""))
  
  bind$feature <- factor(bind$feature, levels = vector$feature)
  bind$meta <- factor(bind$meta, level = paste(listofmeta, "l", sep = ""))
  bind$name<-gsub("[^[:alnum:][:blank:]+?&/\\-]", " ", bind$feature) ## without underbar
  bind$name <- factor(bind$name, levels = bind$name[1:(nrow(bind)/length(listofmeta))])
  
  return (bind)
}

exposures <-  c(
  'bmi_paq1', 'whr_paq1', 'waist_paq1', 'wt_paq', 'wtchg', 
  'logcrp_plasma', 'loghdlc_plasma', 'logtg_plasma', 'logtc_plasma')

### A) exposures ####
bind_exp <- plot_maaslin(exposures, by = 'omega3_ddr')

ggplot(bind_exp, aes(meta, name)) +
  geom_tile(aes(fill = coef), color = "white") +
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


ggplot(data = NULL, aes(species$s__Eubacterium_rectale, totmeta_stn$dpa_ddr))+
  geom_point()


### b) log exposures ####
bind_logexp <- plot_maaslin(logexposures)

ggplot(bind_logexp, aes(meta, name)) +
  geom_tile(aes(fill = coef), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=3, show.legend = TRUE) +
  xlab("log omega3 intake") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "left",
        plot.title = element_text(size=10),
        axis.title.x=element_text(size=10,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10, face="italic",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=10,color="black")) +
  labs(fill = "Beta coefficient")

ggplot(data = NULL, aes(species$s__Eubacterium_rectale, totmeta_stn$dpa_ddr))+
  geom_point()

### c) plasma outcomes ####
bind_out <- plot_maaslin(outcomes, by = 'logcrp_plasma')

ggplot(bind_out, aes(meta, name)) +
  geom_tile(aes(fill = coef), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=3, show.legend = TRUE) +
  xlab("plasma outcomes") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "left",
        plot.title = element_text(size=10),
        axis.title.x=element_text(size=10,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10, face="italic",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=10,color="black")) +
  labs(fill = "Beta coefficient")

ggplot(data = NULL, aes(species$s__Eubacterium_rectale, totmeta_stn$dpa_ddr))+
  geom_point()


### d) exposure with ffq ####
bind_ffq <- plot_maaslin(exposures_ffq, by = 'logcrp_plasma')

ggplot(bind_ffq, aes(meta, name)) +
  geom_tile(aes(fill = coef), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=3, show.legend = TRUE) +
  xlab("plasma outcomes") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "left",
        plot.title = element_text(size=10),
        axis.title.x=element_text(size=10,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10, face="italic",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=10,color="black")) +
  labs(fill = "Beta coefficient")

ggplot(data = NULL, aes(species$s__Eubacterium_rectale, totmeta_stn$dpa_ddr))+
  geom_point()


bind_ffq <- plot_maaslin(c('omega3_noala_ddr', 'epa_ddr', 'dpa_ddr', 'dha_ddr',
                           'omega3_noala_ffq', 'epa_ffq', 'dpa_ffq', 'dha_ffq',
                           'omega3_noala10v', 'epa10v', 'dpa10v', 'dha10v',
                           'logcrp_plasma'), by = 'logcrp_plasma')

bind_ffq <- plot_maaslin(c('logomega3_noala_ddr', 
                           'logomega3_noala_ffq', 
                           'logomega3_noala10v', 
                           'logcrp_plasma'), by = 'logcrp_plasma')

ggplot(bind_ffq, aes(meta, name)) +
  geom_tile(aes(fill = coef), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=3, show.legend = TRUE) +
  xlab("plasma outcomes") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "left",
        plot.title = element_text(size=10),
        axis.title.x=element_text(size=10,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10, face="italic",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=10,color="black")) +
  labs(fill = "Beta coefficient")

ggplot(data = NULL, aes(species$s__Eubacterium_rectale, totmeta_stn$dpa_ddr))+
  geom_point()






### e)  ####

list_maaslin <- c('logomega3_ddr', 'logala_ddr','logepa_ddr' ,'logdpa_ddr' ,'logdha_ddr', 'logomega3_noala_ddr' , 
                  'logala_pfa', 'logepa_pfa', 'logdpa_pfa', 'logdha_pfa', 
                  'logAA_pfa', 'logcrp_plasma', 'loghdl_plasma', 'logtc_plasma', 'logtg_plasma')



bind_list <- plot_maaslin(list_maaslin, by = 'logcrp_plasma')

ggplot(bind_list, aes(meta, name)) +
  geom_tile(aes(fill = coef), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=3, show.legend = TRUE) +
  xlab("") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "left",
        plot.title = element_text(size=10),
        axis.title.x=element_text(size=10,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10, face="italic",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=10,color="black")) +
  labs(fill = "Beta coefficient")

# create right bar for categorizing species to phyla
bind_list$category<-as.factor(bind_list$phylum)
ggplot(bind_list, aes(meta, name))+
  scale_y_discrete(position = "right")+
  geom_tile(aes(fill = category)) +
  xlab("") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title.y= element_blank(),
        axis.text.y = element_blank(),
        axis.title=element_text(size=10,face="bold",colour="white"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=10,colour="white")) +
  labs(fill = "Phyla")








### A) exposures ####
exposures <- c('bmi10', 'wtchg', 'waist_paq1', 'waist'
               'logcrp_plasma', #'omega3_avg', 'omega6_avg',  'trans_avg', 
               'tfat_avg', 
               'omega3_pfa', 'omega6_pfa','sat_pfa', 'monounsat_pfa', 'trans_pfa' )

exposures <- c('bmi10', 'wt_paq', #'height_paq1', 
               'bmi_paq1', 'waist_paq1', #'hip_paq1', 
               'whr_paq1', 'wtchg')
getwd()


bind_exp <- plot_maaslin(exposures, by = 'omega3_ddr')

ggplot(bind_exp, aes(meta, name)) +
  geom_tile(aes(fill = coef), color = "white") +
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


ggplot(data = NULL, aes(species$s__Eubacterium_rectale, totmeta_stn$dpa_ddr))+
  geom_point()













####***adjusted****#


bind_adj <- plot_maaslin(adjust, by = 'adjust')

sig_results_adj<- function(exposure){
  dir = paste('./maaslin_results/adjusted_', exposure, '.pcl/significant_results.tsv',sep='')
  result<-read.table(file = dir, sep = '\t', header = TRUE, check.names=FALSE)
  sig_result <- subset(result, metadata==paste(exposure,'_pfa',sep=''), select =c(feature, metadata))
  return(sig_result)
}

all_results_adj <- function(exposure, label){
  dir = paste('./maaslin_results/adjusted_', exposure, '.pcl/all_results.tsv',sep='')
  result<-read.table(file =dir, sep = '\t', header = TRUE, check.names=FALSE)
  all_result <- subset(result, metadata==paste(exposure,'_pfa',sep=''), select =-c(value, N))
  all_result$meta <- label
  return(all_result)
}


  join <- sig_results_adj(adjust[1])
  for (i in 2:length(adjust)) {
    join<-full_join(join, sig_results_adj(adjust[i]), by="feature")
  }
  sigf_exp<-subset(join, select=feature)
  
  bind <- left_join(sigf_exp, all_results_adj(adjust[1], label = paste(adjust[1], "l", sep = "")), by="feature")
  for (i in 2:length(adjust)) {
    bind <- rbind(bind, left_join(sigf_exp, all_results_adj(adjust[i], label = paste(adjust[i], "l", sep = "")), by="feature"))
  }
  
  bind<-left_join(bind, all_species_name, by="feature")
  
  bind$feature <- substring(bind$feature, 4)
  bind$phylum <- substring(bind$phylum, 4)
  
  bind<-bind[order(bind$meta, bind$phylum, bind$feature),]
  vector <- bind %>%
    arrange(feature) %>% arrange(genus) %>% arrange(family) %>% arrange(order) %>% arrange(class) %>% arrange(phylum) %>% arrange(kingdom) %>% arrange(meta) %>%
    slice(1:(nrow(bind)/length(adjust)))
  
  rownames(bind) <- 1:nrow(bind)
  
  bind$stars <- cut(bind$qval, breaks=c(-Inf, 0.01, 0.05, 0.1, 0.25, Inf), label=c("****", "***", "**", "*", ""))
  
  bind$feature <- factor(bind$feature, levels = vector$feature)
  bind$meta <- factor(bind$meta, level = paste(adjust, "l", sep = ""))
  bind$name<-gsub("[^[:alnum:][:blank:]+?&/\\-]", " ", bind$feature) ## without underbar
  bind$name <- factor(bind$name, levels = bind$name[1:(nrow(bind)/length(adjust))])
  

bind_adj<- bind

ggplot(bind_adj, aes(meta, name)) +
  geom_tile(aes(fill = coef), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=3, show.legend = TRUE) +
  xlab("corr. with plasmic FA, adjusted for FA intake") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "left",
        plot.title = element_text(size=10),
        axis.title.x=element_text(size=10,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10, face="italic",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=10,color="black")) +
  labs(fill = "Beta coefficient")

ggplot(data = NULL, aes(species$s__Eubacterium_rectale, totmeta_stn$dpa_ddr))+
  geom_point()

#####

outcomes <- c('logcrp_plasma', 'loghdl_plasma', 'logtc_plasma', 'logtg_plasma', 
              'AA_pfa', 'ala_pfa', 'epa_pfa', 'dpa_pfa', 'dha_pfa', 'omega3_pfa',
              'logAA_pfa', 'logala_pfa', 'logepa_pfa', 'logdpa_pfa', 'logdha_pfa', 'logomega3_noala_pfa')

pfa <- c('ala_pfa', 'epa_pfa', 'dpa_pfa', 'dha_pfa', 'omega3_pfa','omega3_noala_pfa','omega6_pfa')

bind_pfa <- plot_maaslin(pfa, by = 'logcrp_plasma')

ggplot(bind_pfa, aes(meta, name)) +
  geom_tile(aes(fill = coef), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=3, show.legend = TRUE) +
  xlab("") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "left",
        plot.title = element_text(size=10),
        axis.title.x=element_text(size=10,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10, face="italic",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=10,color="black")) +
  labs(fill = "Beta coefficient")


bindbind <- rbind(bind_pfa, bind_adj)
ggplot(bindbind, aes(meta, name)) +
  geom_tile(aes(fill = coef), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=3, show.legend = TRUE) +
  xlab("") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "left",
        plot.title = element_text(size=10),
        axis.title.x=element_text(size=10,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10, face="italic",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=10,color="black")) +
  labs(fill = "Beta coefficient")
















##### ) graphlan ####
list_graphlan <- c('logomega3_ddr', 'logomega3_noala_ddr' , 
                   'logAA_pfa', 'logcrp_plasma')

bind_graph <- plot_maaslin(list_graphlan, by = 'logcrp_plasma')

ggplot(bind_graph, aes(meta, name)) +
  geom_tile(aes(fill = coef), color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", space = "Lab" ) +
  geom_text(aes(label=stars), color="black", size=3, show.legend = TRUE) +
  xlab("") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "left",
        plot.title = element_text(size=10),
        axis.title.x=element_text(size=10,color="black"),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10, face="italic",color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=10,color="black")) +
  labs(fill = "Beta coefficient")


####### adjusted #####

























##### GraPhlan #####

#+ honestly, we need to color them all... but later.

# Figure 2a: create input data for phylogenetic tree creation GraPhlan >> maybe after.....
tree<-read.delim(file='./data_generated/mlvstree.txt')






bind_graph$coep<-as.integer(bind_graph$coef*1000)

bind_graphcoeppos <- subset(bind_graph, coep> 0)
bind_graphcoepneg <- subset(bind_graph, coep<= 0)
bind_graphcoepposDeciles<-quantile(bind_graphcoeppos$coep, prob = seq(0, 1, length = 11), type = 5)
bind_graphcoepposDeciles
#0%  10%  20%  30%  40%  50%  60%  70%  80%  90% 100% 
#1    1    1    2    3    4    5    6    7    9   25 
bind_graphcoepnegDeciles<-quantile(bind_graphcoepneg$coep, prob = seq(0, 1, length = 11), type = 5)
bind_graphcoepnegDeciles
#  0%  10%  20%  30%  40%  50%  60%  70%  80%  90% 100% 
#-16   -5   -3   -2   -1   -1   -1    0    0    0    0 

bind_graph$coep <- ifelse(bind_graph$coep<(-6), -6, bind_graph$coep)
bind_graph$coep <- ifelse(bind_graph$coep>8, 8, bind_graph$coep)
bind_graph$coep <- ifelse(bind_graph$coep==0 &bind_graph$coef<0, -0.5, bind_graph$coep)
bind_graph$coep <- ifelse(bind_graph$coep==0 &bind_graph$coef>0, 0.5, bind_graph$coep)

#bind_graph<-subset(bind_graph, select=-name)


out_rings_color <-function(meta, num){
  sig_med <- bind_graph[bind_graph$metadata==meta,]
  sig_med$species<-tolower(sig_med$feature)
  match<-full_join(all_species_name, sig_med, by="feature")
  
  match$anno2<-"ring_color"
  match$ring<-num
  match$color[match$coep==0.5]<-"#c1e0d1"
  match$color[match$coep==1]  <-"#b2d8c5"
  match$color[match$coep==2]  <-"#a3d0ba"
  match$color[match$coep==3]  <-"#93c9ae"
  match$color[match$coep==4]  <-"#84c1a3"
  match$color[match$coep==5]  <-"#75b997"
  match$color[match$coep==6]  <-"#66b28c"
  match$color[match$coep==7]  <-"#5ba07e"
  match$color[match$coep==8]  <-"#518e70"
  
  match$color[match$coep==-0.5]<-"#ffc5c5"
  match$color[match$coep==-1]  <-"#ffb2b2"
  match$color[match$coep==-2]  <-"#ff9f9f"
  match$color[match$coep==-3]  <-"#ff8c8c"
  match$color[match$coep==-4]  <-"#ff7979"
  match$color[match$coep==-5]  <-"#ff6666"
  match$color[match$coep==-6]  <-"#ff5353"
  
  match$color[is.na(match$coep)]<-"#CCCCC3"
  match$color<-toupper(match$color)
  
  color<-subset(match, select = c(name, anno2, ring, color))
  
  return(color)
}

co1<-out_rings_color(meta="logomega3_ddr", num=7)
co2<-out_rings_color(meta = "logomega3_noala_ddr" , num=6)
co3<-out_rings_color(meta = "logAA_pfa", num=5)
co4<-out_rings_color(meta = "logcrp_plasma", num=4)


coall <-rbind(co1, co2, co3, co4)
write.table(coall, "./data/out_rings_color.txt", 
            sep="\t", quote = FALSE, col.names = F, row.names = F)

out_rings_alpha <-function(path, meta, num){
  
  sig_med <- bind_graph[bind_graph$metadata==meta,]
  
  sig_med$feature<-tolower(sig_med$feature)
  
  match<-full_join(all_species_name, sig_med, by="feature")
  
  match$anno1<-"ring_alpha"
  match$ring<-num
  match$alpha[is.na(match$coef)]<-0.3
  match$alpha[!is.na(match$coef)]<-0.8
  
  alpha<-subset(match, select=c(name, anno1, ring, alpha))
  return(alpha)
}


al1<-out_rings_alpha(path="./maaslin_results/logomega3_ddr.pcl/significant_results.tsv", meta="paee_pam", num=7)
al2<-out_rings_alpha(path="./maaslin_results/logomega3_noala_ddr.pcl/significant_results.tsv",    meta = "act_paqlong", num=6)
al3<-out_rings_alpha(path="./maaslin_results/logAA_pfa.pcl/significant_results.tsv",   meta = "bmi_dlw", num=5)
al4<-out_rings_alpha(path="./maaslin_results/logcrp_plasma.pcl/significant_results.tsv",     meta = "pfat_dlw", num=4)


alall<-rbind(al1, al2, al3, al4)
write.table(alall, "./data/out_rings_alpha.txt", 
            sep="\t", quote = FALSE, col.names = F, row.names = F)
write.table(all_species_name, "./data/all_species_name.txt", 
            sep=".", quote = FALSE, col.names = F, row.names = F)
