
##############################################################################################################################################
# 1) Purposes: conduct PERMANOVA using brey-curtis to calculation the contribution from the exposure, outcome, and covariates to the variarion of overall microbial structure
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

library(ggplot2)
library(scales)
library(tidyverse)
library(vegan)

### IMPORT DATA #####
##### + with NAs ####
species <-read.table(file = './data/species_filt.csv',
                     sep = ',',  row.names = 1,  header = TRUE,   check.names = FALSE)

meta_stn <- read.table('./data/meta_stn.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")

totmeta_stn <- read.table('./data/totmeta_stn.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")


# calculate bray-curtis matrix
## don't think it is needed in adonis2, was needed in adonis (1)
# tspecies_bc <- vegdist(t(species), "bray")

# PERMANOVA using adonis2, strata=as.factor(meta_med$SubjectID)
colnames(meta_stn)
meta_stn$ala_pfa
perm_list <- list(omega3_w, omega3_noala_w, ala_w, epa_w, dha_w, dpa_w,
                  omega3_avg, omega3_noala_avg, ala_avg, epa_avg, dpa_avg, dha_avg,
                  ala_pfa, epa_pfa, dha_pfa, dpa_pfa, # only 1 plasma metabolite data
                  bmi10, calor10n, wt10, act10,calor_fs_dr_w,#fatmass, 6moweightchange, changesinceage21,
                  agemlvs, smoke10, alco10m, abx_avg, alc_avg, probx_avg, bristol_avg,
                  
                  logcrp_w, crp_w, crp_avg,
                  hdlc_w, tchdl_w, tc_w, 
                  logtchdl_w, logtc_w, loghdl_w, logtg_w,tg_w 
                  )
perm_list <- c("omega3_w", "omega3_noala_w", "ala_w", "epa_w", "dha_w", "dpa_w",
                  'omega3_avg', 'omega3_noala_avg', 'ala_avg', 'epa_avg', "dpa_avg", "dha_avg",
                  'ala_pfa', 'epa_pfa', 'dha_pfa', 'dpa_pfa', # only 1 plasma metabolite data
                  'bmi10', 'calor10n', 'wt10', 'act10','calor_fs_dr_w',#fatmass, 6moweightchange, changesinceage21,
                  'agemlvs', 'smoke10','alco10n', 'abx_avg', 'alc_avg', 'probx_avg', 'bristol_avg',
                  'logcrp_w', 'crp_w', 'crp_avg',
                  'hdlc_w', 'tchdl_w', 'tc_w', 
                  'logtchdl_w', 'logtc_w', 'loghdl_w', 'logtg_w','tg_w' )
length(perm_list)
#permanova <- data.frame()

# meta_stn[,perm_list[1]]
p=perm_list[2]
# as.name(p)


permanova<- NULL


meta_stn_num <- colnames(meta_stn[ , unlist(lapply(meta_stn, is.numeric)) ])



for (p in meta_stn_num){
   s <- species
   m <- meta_stn
   a <- adonis2(formula = s ~ m[,p], strata=m$ID1,  
                data = m, permutations = 999, na.action = na.exclude)
   x <- data.frame(name = p, a[1,c("R2","Pr(>F)")])
   permanova <- rbind (permanova, x)
}

row.name(permanova) <- perm_list


#+ na.action = 
#+ na.fail ; defalt
#+ na.omit ; remove cases
#+ na.exclude ; differs from omit in na.action, used -> padded to correct length by inserting NAs

#+ Trouble shooting
#+ 1) ['qr' and 'y' must have the same number of rows]
#+ species ~ logcrp_w don't have the same number of rows
#+ 2) [Number of observations and length of Block 'strata' do not match.]
#+ there shouldn't be any NAs in logcrp_w


permanova$name <- factor(permanova$name, levels = rev(colnames(meta_stn[ , unlist(lapply(meta_stn, is.numeric)) ])))
permanova$R2 <- as.numeric(permanova$R2)
permanova$stars <- cut(permanova$Pr..F., breaks=c(-Inf, 0.01, 0.05, 0.10, Inf), label=c("***","**", "*", ""))
permanova <- permanova %>%
  mutate(colors = ifelse(name %in% basic, "basic", 
                         ifelse()))

write.table(permanova, file = "data/permanova.tsv", sep = "\t", col.names = TRUE, row.names = TRUE)



## list comes from Fig0
basic <-  data.frame(name = c('ID1_stn', 'ID1', 'phase',
             'agemlvs', 'bmi10', 'bmi12', 'wt10', 'wt12', 'smoke10', 'smk12', 'ncig12', 'bmic10', 'bmic12', 'bmi12cat', 'smoke12', 
             'calor10n', 'alco10n', 'calor10v', 'calor_avg', 'alco_avg', 'abx_avg', 'probx_avg', 'alc_avg', 'bristol', 'bristolcat', 'bristol_avg', 'bristolcat_avg',
             'a_alco_fo_dr_ddr', 'abx_ddr', 'probx_ddr', 'bristol_ddr', 'bristolcat_ddr' ), type = "basic")
# smoke12 is ordinal
exposure_long <- data.frame(name =c('ala10v', 'epa10v', 'dha10v', 'trans10v', 'omega610v', 'omega3_noala10a', 
                                    'dpa10v', 'omega310v', 'omega3_noala10v'), type = "exposure_long")

exposure_avg <- data.frame(name =c('ala_avg', 'epa_avg', 'dpa_avg', 'trans_avg', 
                                   'omega3_avg', 'omega6_avg', 'omega3_noala_avg', 'fishavg'), type = "exposure_avg")

exposure_short <- data.frame(name =c( 'ala_ddr',
                     'epa_ddr' ,
                     'dha_ddr',
                     'dpa_ddr' ,
                     'omega3_noala_ddr',
                     'omega3_ddr',
                     'a_fat_fs_dr_ddr'), type = "exposure_short")

exposure_1yr <- data.frame(name =c('ala_ffq', 'epa_ffq', 'dpa_ffq', 'dha_ffq', 'omega3_ffq', 'omega6_ffq', 'omega3_noala_ffq', 'omega_3_epa_ffq',
                  'tuna_ffq', 'dk_fish_ffq', 'oth_fish_ffq','fishffq'), type = "exposure_1yr")

outcome_plasma <- data.frame(name =c('folate_plasma', 'tc_plasma', 'hdlc_plasma', 'tg_plasma', 'crp_plasma',  'crp_plasmacat', 'tchdl_plasma', 
                    'logtchdl_plasma', 'logcrp_plasma', 'logtc_plasma', 'loghdl_plasma', 'logtg_plasma'), type = "outcome_plasma")

outcome_pfa <- data.frame(name =c('ala_pfa', 'epa_pfa', 'dpa_pfa', 'dha_pfa', 
                 'AA_pfa' , 
                 'omega3_pfa','omega3_noala_pfa','omega6_pfa'), type = "outcome_pfa")
types <- rbind(basic, exposure_long, exposure_avg, exposure_short, exposure_1yr, outcome_plasma, outcome_pfa)

permanova2 <- permanova %>% filter(!(name %in% c("ID1", "omega_3_epa_ffq")))
permanova2 <- left_join(permanova2, types, by = "name")
permanova2$name <- factor(permanova2$name, levels = rev(colnames(meta_stn[ , unlist(lapply(meta_stn, is.numeric)) ])),
                          ordered = T)

ggplot(permanova2, aes(x=R2, y=name, fill = type)) +
  geom_bar(stat="identity")+
  theme_light()+
  scale_x_continuous(labels=percent, limits = c(0, 0.01))+
  geom_text(aes(label=stars), color="black", size=5) +
  theme_classic()+
  theme(axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"),
        legend.position = "none",
        legend.text = element_text(size = 10,color="black"),
        legend.title = element_blank(),
        plot.title = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10,color="black"),
        axis.text.x = element_text(size=10,color="black"))  


histogram(log(meta_stn$ala_pfa))





#### LOG PERMANOVA #####

logpermanova<- NULL
logmeta_stn <- read.table('./data/logmeta_stn.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")


logmeta_stn_num <- colnames(logmeta_stn[ , unlist(lapply(logmeta_stn, is.numeric)) ])

for (p in logmeta_stn_num){
  s <- species
  m <- logmeta_stn
  a <- adonis2(formula = s ~ m[,p], strata=m$ID1,  
               data = m, permutations = 999, na.action = na.exclude)
  x <- data.frame(name = p, a[1,c("R2","Pr(>F)")])
  logpermanova <- rbind (logpermanova, x)
}



logpermanova$R2 <- as.numeric(logpermanova$R2)
logpermanova$stars <- cut(logpermanova$Pr..F., breaks=c(-Inf, 0.01, 0.05, 0.10, 0.20, Inf), label=c("***","**", "*", ".", ""))


write.table(logpermanova, file = "data/logpermanova.tsv", sep = "\t", col.names = TRUE, row.names = TRUE)



## list comes from Fig0
logtypes <- types
logtypes$name <- paste("log", types$name, sep = "")

# logpermanova2 <- logpermanova %>% filter(!(name %in% c("ID1", "omega_3_epa_ffq")))
logpermanova2 <- left_join(logpermanova, logtypes, by = "name")
logpermanova2$name <- factor(logpermanova2$name, levels = rev(colnames(logmeta_stn[ , unlist(lapply(logmeta_stn, is.numeric)) ])),
                          ordered = T)

ggplot(logpermanova2, aes(x=R2, y=name, fill = type)) +
  geom_bar(stat="identity")+
  theme_light()+
  scale_x_continuous(labels=percent, limits = c(0, 0.01))+
  geom_text(aes(label=stars), color="black", size=5) +
  theme_classic()+
  theme(axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"),
        legend.position = "none",
        legend.text = element_text(size = 10,color="black"),
        legend.title = element_blank(),
        plot.title = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=10,color="black"),
        axis.text.x = element_text(size=10,color="black"))  


histogram(log(meta_stn$epa_ddr))























short= c("omega3_w", "omega3_noala_w", "ala_w", "epa_w", "dha_w", "dpa_w")
long= c('omega3_avg', 'omega3_noala_avg', 'ala_avg', 'epa_avg', "dpa_avg", "dha_avg")
plasma = c('ala_pfa', 'epa_pfa', 'dha_pfa', 'dpa_pfa')
health = c('bmi10', 'calor10n', 'wt10', 'act10','calor_fs_dr_w',#fatmass, 6moweightchange, changesinceage21,
           'agemlvs', 'smoke10','alco10n', 'abx_avg', 'alc_avg', 'probx_avg', 'bristol_avg')
inflammation =  c('logcrp_w', 'crp_w', 'crp_avg',
                  'hdlc_w', 'tchdl_w', 'tc_w', 
                  'logtchdl_w', 'logtc_w', 'loghdl_w', 'logtg_w','tg_w' )





permanova %>% 
  mutate(cat = ifelse(sum(name == short), "short",
                      ifelse(sum(name == long), "long",
                             ifelse(sum(name == plasma), "plasma",
                                    ifelse(sum(name == health), "health","inflammation")))))
                        

all_taxonomy<-cbind(all_taxonomy, component)
all_taxonomy_pa<-all_taxonomy[1:8,]
all_taxonomy_pa$cat<-"PA"

all_taxonomy_wt<-all_taxonomy[9:12,]
all_taxonomy_wt$cat<-"Body weight"

all_taxonomy_bio<-all_taxonomy[13:14,]
all_taxonomy_bio$cat<-"Plasma biomarkers"

all_taxonomy_cov<-all_taxonomy[15:21,]
all_taxonomy_cov$cat<-"Covariables"

all_taxonomy_tax<-rbind(all_taxonomy_pa, all_taxonomy_wt, all_taxonomy_bio, all_taxonomy_cov)
all_taxonomy_tax$cat<-as.factor(all_taxonomy_tax$cat)

color_code<-c("PA"="#3EAB07",
              "Body weight"="#FF1713",
              "Plasma biomarkers"="#FFB713",
              "Covariables"="#8dd3c7")

all_taxonomy_tax$stars <- cut(all_taxonomy_tax$`Pr(>F)`, breaks=c(-Inf, 0.01, 0.05, 0.10, Inf), label=c("***","**", "*", ""))
write.csv(all_taxonomy_tax, file="./data_generated/permanova_tax.csv")

all_taxonomy_tax<-read.csv(file="./data_generated/permanova_tax.csv")

level_y_order <- factor(all_taxonomy_tax$component, level = rev(all_taxonomy_tax$component))

# create figure to show PERMANOVA results
png(file="./figure_generated/fig2b_barplot_permanova.png",width=2600, height=2000, pointsize=50)
ggplot(all_taxonomy_tax, aes(level_y_order, R2, fill=cat)) +
  geom_bar(stat="identity")+
  theme_light()+scale_y_continuous(labels=percent, limits = c(0, 0.01))+
  geom_text(aes(label=stars), color="black", size=40) +
  scale_fill_manual(values = color_code)+ coord_flip()+
  theme(axis.line = element_line(colour = "black", 
                                 size = 2, linetype = "solid"),
        legend.position = "none",
        legend.text = element_text(size = 60,color="black"),
        legend.title = element_blank(),
        plot.title = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y= element_blank(),
        axis.text.y=element_text(size=60,color="black"),
        axis.text.x = element_text(size=60,color="black"))  
dev.off()
