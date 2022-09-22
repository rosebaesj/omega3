
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
setwd("/udd/nhkwa/mlvs")
options(max.print=1000000)

library(Maaslin2)
library(ggplot2)
library(grid)
library(ggpubr)
library(tidyverse)
library(gridExtra)
library(colourvalues)

# read in metadata
meta_med<-read.csv(file="./data_generated/meta_species.csv")
meta_med$mvpa_paqlong <- meta_med$vpa_paqlong + meta_med$mpa_paqlong

summary(meta_med$act_paqlong_sec)
summary(meta_med$act_paqlong_como)

# read in taxonomy data
tax_rpk_name <-   read.table(file = './data_generated/bugs_dna_929_unFilt.tsv',
                             sep = '\t',    header = TRUE,    check.names = FALSE,    na.strings = c("", "NA"))
tax_rpk_name<-tax_rpk_name %>%
  separate(Sample, c("kingdom",       "phylum",        "class" ,        "order",         "family",        "genus" ,        "species" ,      "strain"), 
           sep = '\\|', remove = TRUE)


top20 <- colnames(species_filt[rev(order(colSums(species_filt)))[1:20]])
meta_stn$omega3_w

meta_species_nona$logcrp_w

ggplot(meta_species_nona, aes(s__Bacteroides_uniformis, log(omega3_noala_w))) + 
  geom_point(aes(color = logcrp_w),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff",  "#dce319ff","#b8de29ff","#29af7fff", "#2d708eff", "#440154ff"))+
  #xlab(expression(paste('Relative abundance of',italic("A. putredinis"), '(log10 scale)')))+
  #ylab("Recent total PA (MET-hours/week)")+
  theme_classic()+
  theme(legend.position="right",
        axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"),
        legend.text = element_text(size=10),
        #legend.title = element_blank(),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10),
        axis.title.y=element_text(color = "Black",size=10),
        #axis.title.x=element_blank()
        )

top20



logCRP ~ omega_noala + BUG +omega_noala :BUG 
+age +antibiotic +calorie + (1|participant)


lm(data = )

meta_species_nona$abx_w
summary(lm(logcrp_w ~ log(omega3_noala_w) + s__Bacteroides_uniformis
             agemlvs +abx +calor_fs_dr +(1|participant), data = meta_species_nona))


top20

summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Bacteroides_uniformis +log(omega3_noala_w)*s__Bacteroides_uniformis 
             +agemlvs + abx_w+(1|ID1)))
sum(meta_species_nona$s__Bacteroides_uniformis == 0)


summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Bacteroides_vulgatus +log(omega3_noala_w)*s__Bacteroides_vulgatus 
             +agemlvs+ calor_fs_dr_w+ abx_w+(1|ID1)))
sum(meta_species_nona$s__Bacteroides_vulgatus == 0)



summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Faecalibacterium_prausnitzii +log(omega3_noala_w)*s__Faecalibacterium_prausnitzii 
             +agemlvs+ calor_fs_dr_w+ abx_w+(1|ID1)))
sum(meta_species_nona$s__Faecalibacterium_prausnitzii == 0)




summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Eubacterium_rectale +log(omega3_noala_w)*s__Eubacterium_rectale 
             +agemlvs+ calor_fs_dr_w+ abx_w+(1|ID1)))

summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Prevotella_copri +log(omega3_noala_w)*s__Prevotella_copri 
             +agemlvs+ calor_fs_dr_w+ abx_w+(1|ID1)))
sum(meta_species_nona$s__Prevotella_copri == 0)




#5
summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Alistipes_putredinis +log(omega3_noala_w)*s__Alistipes_putredinis 
             +agemlvs+ calor_fs_dr_w+ abx_w+(1|ID1)))

sum(meta_species_nona$s__Alistipes_putredinis == 0)

summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Bacteroides_dorei +log(omega3_noala_w)*s__Bacteroides_dorei 
             +agemlvs+ calor_fs_dr_w+ abx_w+(1|ID1)))
sum(meta_species_nona$s__Bacteroides_dorei == 0)


summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Eubacterium_eligens +log(omega3_noala_w)*s__Eubacterium_eligens 
             +agemlvs+ calor_fs_dr_w+ abx_w+(1|ID1)))
sum(meta_species_nona$s__Eubacterium_eligens == 0)

summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Roseburia_faecis +log(omega3_noala_w)*s__Roseburia_faecis 
             +agemlvs+ calor_fs_dr_w+ abx_w+(1|ID1)))

summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Bacteroides_stercoris +log(omega3_noala_w)*s__Bacteroides_stercoris 
             +agemlvs+ calor_fs_dr_w+ abx_w+(1|ID1)))

ggplot(meta_species_nona, aes(s__Eubacterium_siraeum, log(omega3_noala_w))) + 
  geom_point(aes(color = logcrp_w),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff",  "#dce319ff","#b8de29ff","#29af7fff", "#2d708eff", "#440154ff"))+
  #xlab(expression(paste('Relative abundance of',italic("A. putredinis"), '(log10 scale)')))+
  #ylab("Recent total PA (MET-hours/week)")+
  theme_classic()+
  theme(legend.position="right",
        axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"),
        legend.text = element_text(size=10),
        #legend.title = element_blank(),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10),
        axis.title.y=element_text(color = "Black",size=10),
        #axis.title.x=element_blank()
  )
summary(meta_species_nona$s__Eubacterium_siraeum)


sum(meta_species_nona$s__Eubacterium_siraeum != 0)

##11


summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Eubacterium_siraeum +log(omega3_noala_w)*s__Eubacterium_siraeum 
             +agemlvs+ calor_fs_dr_w+ abx_w+(1|ID1)))
sum(meta_species_nona$s__Eubacterium_siraeum == 0)
ggplot(meta_species_nona, aes(s__Eubacterium_siraeum, log(omega3_noala_w))) + 
  geom_point(aes(color = logcrp_w),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff",  "#dce319ff","#b8de29ff","#29af7fff", "#2d708eff", "#440154ff"))+
  #xlab(expression(paste('Relative abundance of',italic("A. putredinis"), '(log10 scale)')))+
  #ylab("Recent total PA (MET-hours/week)")+
  theme_classic()+
  theme(legend.position="right",
        axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"),
        legend.text = element_text(size=10),
        #legend.title = element_blank(),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10),
        axis.title.y=element_text(color = "Black",size=10),
        #axis.title.x=element_blank()
  )

summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Ruminococcus_bromii +log(omega3_noala_w)*s__Ruminococcus_bromii 
             +agemlvs+ calor_fs_dr_w+ abx_w+(1|ID1)))
sum(meta_species_nona$s__Ruminococcus_bromii == 0)
ggplot(meta_species_nona, aes(s__Ruminococcus_bromii, log(omega3_noala_w))) + 
  geom_point(aes(color = logcrp_w),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff",  "#dce319ff","#b8de29ff","#29af7fff", "#2d708eff", "#440154ff"))+
  #xlab(expression(paste('Relative abundance of',italic("A. putredinis"), '(log10 scale)')))+
  #ylab("Recent total PA (MET-hours/week)")+
  theme_classic()+
  theme(legend.position="right",
        axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"),
        legend.text = element_text(size=10),
        #legend.title = element_blank(),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10),
        axis.title.y=element_text(color = "Black",size=10),
        #axis.title.x=element_blank()
  )

summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Fusicatenibacter_saccharivorans +log(omega3_noala_w)*s__Fusicatenibacter_saccharivorans 
             +agemlvs+ calor_fs_dr_w+ abx_w+(1|ID1)))

summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Parabacteroides_distasonis +log(omega3_noala_w)*s__Parabacteroides_distasonis 
             +agemlvs+ calor_fs_dr_w+ abx_w+(1|ID1)))

summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Akkermansia_muciniphila +log(omega3_noala_w)*s__Akkermansia_muciniphila 
             +agemlvs+ calor_fs_dr_w+ abx_w+(1|ID1)))

summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Bacteroides_eggerthii +log(omega3_noala_w)*s__Bacteroides_eggerthii 
             +agemlvs+ calor_fs_dr_w+ abx_w+(1|ID1)))
sum(meta_species_nona$s__Bacteroides_eggerthii == 0)
ggplot(meta_species_nona, aes(s__Bacteroides_eggerthii, log(omega3_noala_w))) + 
  geom_point(aes(color = logcrp_w),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff",  "#dce319ff","#b8de29ff","#29af7fff", "#2d708eff", "#440154ff"))+
  #xlab(expression(paste('Relative abundance of',italic("A. putredinis"), '(log10 scale)')))+
  #ylab("Recent total PA (MET-hours/week)")+
  theme_classic()+
  theme(legend.position="right",
        axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"),
        legend.text = element_text(size=10),
        #legend.title = element_blank(),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10),
        axis.title.y=element_text(color = "Black",size=10),
        #axis.title.x=element_blank()
  )




summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Ruminococcus_bicirculans +log(omega3_noala_w)*s__Ruminococcus_bicirculans 
             +agemlvs+ calor_fs_dr_w+ abx_w+(1|ID1)))

summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Bacteroides_cellulosilyticus +log(omega3_noala_w)*s__Bacteroides_cellulosilyticus 
             +agemlvs+ calor_fs_dr_w+ abx_w+(1|ID1)))
summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Bacteroides_cellulosilyticus +log(omega3_noala_w)*s__Bacteroides_cellulosilyticus 
             +agemlvs+ calor_fs_dr_w+ abx_w+(1|ID1)))

sum(meta_species_nona$s__Bacteroides_cellulosilyticus == 0)
ggplot(meta_species_nona, aes(s__Bacteroides_cellulosilyticus, log(omega3_noala_w))) + 
  geom_point(aes(color = logcrp_w),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff",  "#dce319ff","#b8de29ff","#29af7fff", "#2d708eff", "#440154ff"))+
  #xlab(expression(paste('Relative abundance of',italic("A. putredinis"), '(log10 scale)')))+
  #ylab("Recent total PA (MET-hours/week)")+
  theme_classic()+
  theme(legend.position="right",
        axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"),
        legend.text = element_text(size=10),
        #legend.title = element_blank(),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10),
        axis.title.y=element_text(color = "Black",size=10),
        #axis.title.x=element_blank()
  )




summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Bacteroides_thetaiotaomicron +log(omega3_noala_w)*s__Bacteroides_thetaiotaomicron 
             +agemlvs+ calor_fs_dr_w+ abx_w+(1|ID1)))

summary(lmer(data = meta_species_nona, 
             logcrp_w ~ log(omega3_noala_w) + s__Anaerostipes_hadrus +log(omega3_noala_w)*s__Anaerostipes_hadrus 
             +agemlvs+ calor_fs_dr_w+ abx_w+(1|ID1)))


meta_species$smoke12




top20


ggplot(meta_species_nona, aes(s__Bacteroides_uniformis, epa)) + 
  geom_point(aes(color = log(omega3_w)),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff",  "#dce319ff","#b8de29ff","#29af7fff", "#2d708eff", "#440154ff"))+
  #xlab(expression(paste('Relative abundance of',italic("A. putredinis"), '(log10 scale)')))+
  #ylab("Recent total PA (MET-hours/week)")+
  theme_classic()+
  theme(legend.position="right",
        axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10),
        axis.title.y=element_text(color = "Black",size=10),
        #axis.title.x=element_blank()
  )

# scatterplot with x-axis of s__Alistipes_putredinis and y-axis of tdee_pam colored by bmi_dlw
png(file="./figure_generated/pa_bmi_Alistipes.png",width=2000,height=2000, pointsize=80)
ggplot(ttax_rel, aes(s__Alistipes_putredinis.y, paee_pamwk)) + geom_point(aes(color = bmi_dlw),size=20, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff",  "#dce319ff","#b8de29ff","#29af7fff", "#2d708eff", "#440154ff"))+
  xlab(expression(paste('Relative abundance of',italic("A. putredinis"), '(log10 scale)')))+
  ylab("Recent total PA (MET-hours/week)")+theme_classic()+
  theme(legend.position="right",
        axis.line = element_line(colour = "black", 
                                 size = 2, linetype = "solid"),
        legend.text = element_text(size=80),
        legend.title = element_blank(),
        axis.text.y=element_text(color = "Black",size=80),
        axis.text.x=element_text(color = "Black",size=80),
        axis.title.y=element_text(color = "Black",size=80),
        axis.title.x=element_blank())
dev.off()

# scatterplot with x-axis of s__Alistipes_putredinis and y-axis of tdee_pam colored by pfat_dlw
png(file="./figure_generated/pa_pfat_Alistipes.png",width=2000,height=2000, pointsize=80)
ggplot(ttax_rel, aes(s__Alistipes_putredinis.y, paee_pamwk)) + geom_point(aes(color = pfat_dlw),size=20, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff",  "#dce319ff","#b8de29ff","#29af7fff", "#2d708eff", "#440154ff"))+
  xlab("Abundance of Alistipes putredinis")+
  ylab("Recent total PA (MET-hours/week)")+theme_classic()+
  theme(legend.position="right",
        axis.line = element_line(colour = "black", 
                                 size = 2, linetype = "solid"),
        legend.text = element_text(size=80),
        legend.title = element_blank(),
        axis.text.y=element_text(color = "Black",size=80),
        axis.text.x=element_text(color = "Black",size=80),
        axis.title.y=element_text(color = "Black",size=80),
        axis.title.x=element_text(color = "Black",size=80))
dev.off()

# scatterplot with x-axis of s__Alistipes_putredinis and y-axis of act_paqlong colored by wtchgsto21
#"#fde725ff",  "#dce319ff","#b8de29ff","#29af7fff", "#2d708eff", "#440154ff"
#"#98ff98","#7be27d","#5fc663","#41aa4a","#005b00","#004200","#002d00","#000d00"
#colour_values(1:10, n_summaries = 6)
#df <- data.frame(a = 10, x = 1:10)
#df$col <- colour_values(-df$x, palette = "viridis")
png(file="./figure_generated/pa_sec_wtchgsto21_Alistipes.png",width=2000,height=2000, pointsize=80)
ggplot(ttax_rel, aes(s__Alistipes_putredinis.y, act_paqlong_sec)) + geom_point(aes(color = wtchgsto21),size=20, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#FDE725FF", "#B4DD2CFF","#5DC963FF","#21908CFF","#3B528BFF","#482878FF","#440154FF"))+
  ylim(0, 90)+
  xlab(expression(paste('Relative abundance of',italic("A. putredinis"), '(log10 scale)')))+
  ylab("Long-term PA (MET-hours/week)")+theme_classic()+
  theme(legend.position="right",
        axis.line = element_line(colour = "black", 
                                 size = 2, linetype = "solid"),
        legend.text = element_text(size=80),
        legend.title = element_blank(),
        axis.text.y=element_text(color = "Black",size=80),
        axis.text.x=element_text(color = "Black",size=80),
        axis.title.y=element_text(color = "Black",size=80),
        axis.title.x=element_blank())
dev.off()

