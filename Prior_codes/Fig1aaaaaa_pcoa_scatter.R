
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
library(tidyverse)
library(vegan)
library(RColorBrewer)
library(grid)
library(ggpubr)
library(viridis)

# read in meta data
totmeta_stn <- read.table('./data/totmeta_stn.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")

# read in taxonomy data
# only keep species-level features

species_filt <- read.table(file = "species_filt.tsv", sep = "\t",header=TRUE, check.names=FALSE)


rowSums(species_filt) 
#+ near 1 but not 1, because calculated relative abundance before filtering. 
#+ this is different from Kai's approach 

# calculate bray-crutis matrix
species_filt_bc <- vegdist(species_filt, "bray")

# PCOA using cmdscale
mod <- cmdscale(species_filt_bc, eig = T)
pcoap <- data.frame(mod$points)
pcoap$id_st<-rownames(pcoap)

# percentages of variation explained by PCO1 & 2
mod$eig[1]/sum(mod$eig)
mod$eig[2]/sum(mod$eig)


meta_pcoa <- merge(pcoap, totmeta_stn, by = 0)









# scatterplot with x-axis of s__Alistipes_putredinis and y-axis of tdee_pam colored by bmi_dlw


# omega3 intake -> log crp

histogram(meta_pcoa$X1)
reg <- lmer( logcrp_plasma ~ logomega3_noala_ddr + X1 + logomega3_noala_ddr*X1 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_pcoa)
summary(reg)
a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
anova(reg)[nrow(anova(reg)),ncol(anova(reg))]


ggplot(meta_pcoa, aes(X1, logomega3_noala_ddr)) + 
  geom_point(aes(color = logcrp_plasma),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff", "#b8de29ff","#29af7fff", "#2d708eff", "#440154ff", "#440154ff"))+
  labs( title = "color by logcrp_plasma",
    subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
  xlab("PCo 1")+
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


a<-quantile(meta_pcoa$logcrp_plasma, probs = seq(0, 1, 0.2))
a[1] <- -Inf
a[6] <- +Inf
meta_pcoa$logcrp_plasmaq <- cut(meta_pcoa$logcrp_plasma, breaks=a, label=c(1, 2, 3, 4, 5))

ggplot(meta_pcoa, aes(x = logcrp_plasmaq, y = logomega3_noala_ddr, fill=logcrp_plasmaq)) +
  geom_boxplot(colour = "black",lwd=0.5, outlier.color = "black",outlier.size = 1)+
  theme_bw() +
 # scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8", "#253494"))+
  ylab("logomega3_noala_ddr")+
  xlab("logcrp_plasmaq")+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=1),
        plot.title=element_blank(),
        axis.title.y=element_text(color = "Black",size=10),
        axis.title.x=element_text(color = "Black",size=10),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10,angle = 300,hjust = 0,vjust=1))



#omega3 intake -> omega3 pfa


reg <- lmer( logomega3_noala_pfa ~ logomega3_noala_ddr + X1 + logomega3_noala_ddr*X1 
             #+agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_pcoa)
summary(reg)
a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
anova(reg)[nrow(anova(reg)),ncol(anova(reg))]

ggplot(meta_pcoa, aes(X1, logomega3_noala_ddr)) + 
  geom_point(aes(color = logomega3_noala_pfa),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff", "#b8de29ff","#29af7fff", "#2d708eff", "#440154ff", "#440154ff"))+
  labs( title = "color by logomega3_noala_pfa",
        subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
  #xlab(top20[i])+
  #ylab("Recent total PA (MET-hours/week)")+
  xlab("PCo 1")+
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

a<-quantile(meta_pcoa$logomega3_noala_pfa, probs = seq(0, 1, 0.2))
a[1] <- -Inf
a[6] <- +Inf
meta_pcoa$logomega3_noala_pfaq <- cut(meta_pcoa$logomega3_noala_pfa, breaks=a, label=c(1, 2, 3, 4, 5))

ggplot(meta_pcoa, aes(x = logomega3_noala_pfaq, y = logomega3_noala_ddr, fill=logomega3_noala_pfaq)) +
  geom_boxplot(colour = "black",lwd=0.5, outlier.color = "black",outlier.size = 1)+
  theme_bw() +
  # scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8", "#253494"))+
  ylab("logomega3_noala_ddr")+
  xlab("logomega3_noala_pfaq")+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=1),
        plot.title=element_blank(),
        axis.title.y=element_text(color = "Black",size=10),
        axis.title.x=element_text(color = "Black",size=10),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10,angle = 300,hjust = 0,vjust=1))


# omega3 intake -> log crp


reg <- lmer( logcrp_plasma ~ omega3_noala_ddr + X1 + omega3_noala_ddr*X1 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_pcoa)
summary(reg)
a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
anova(reg)[nrow(anova(reg)),ncol(anova(reg))]


ggplot(meta_pcoa, aes(X1, omega3_noala_ddr)) + 
  geom_point(aes(color = logcrp_plasma),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff", "#b8de29ff","#29af7fff", "#2d708eff", "#440154ff", "#440154ff"))+
  labs( title = "color by logcrp_plasma",
        subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
  xlab("PCo 1")+
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


a<-quantile(meta_pcoa$logcrp_plasma, probs = seq(0, 1, 0.2))
a[1] <- -Inf
a[6] <- +Inf
meta_pcoa$logcrp_plasmaq <- cut(meta_pcoa$logcrp_plasma, breaks=a, label=c(1, 2, 3, 4, 5))

ggplot(meta_pcoa, aes(x = logcrp_plasmaq, y = omega3_noala_ddr, fill=logcrp_plasmaq)) +
  geom_boxplot(colour = "black",lwd=0.5, outlier.color = "black",outlier.size = 1)+
  theme_bw() +
  # scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8", "#253494"))+
  ylab("omega3_noala_ddr")+
  xlab("logcrp_plasmaq")+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=1),
        plot.title=element_blank(),
        axis.title.y=element_text(color = "Black",size=10),
        axis.title.x=element_text(color = "Black",size=10),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10,angle = 300,hjust = 0,vjust=1))



# omega3 with ala intake -> log crp


reg <- lmer( logcrp_plasma ~ logomega3_ddr + X1 + logomega3_ddr*X1 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_pcoa)
summary(reg)
a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
anova(reg)[nrow(anova(reg)),ncol(anova(reg))]


ggplot(meta_pcoa, aes(X1, logomega3_ddr)) + 
  geom_point(aes(color = logcrp_plasma),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff", "#b8de29ff","#29af7fff", "#2d708eff", "#440154ff", "#440154ff"))+
  labs( title = "color by logcrp_plasma",
        subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
  xlab("PCo 1")+
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


a<-quantile(meta_pcoa$logcrp_plasma, probs = seq(0, 1, 0.2))
a[1] <- -Inf
a[6] <- +Inf
meta_pcoa$logcrp_plasmaq <- cut(meta_pcoa$logcrp_plasma, breaks=a, label=c(1, 2, 3, 4, 5))

ggplot(meta_pcoa, aes(x = logcrp_plasmaq, y = logomega3_ddr, fill=logcrp_plasmaq)) +
  geom_boxplot(colour = "black",lwd=0.5, outlier.color = "black",outlier.size = 1)+
  theme_bw() +
  # scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8", "#253494"))+
  ylab("omega3_noala_ddr")+
  xlab("logcrp_plasmaq")+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=1),
        plot.title=element_blank(),
        axis.title.y=element_text(color = "Black",size=10),
        axis.title.x=element_text(color = "Black",size=10),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10,angle = 300,hjust = 0,vjust=1))

# omega3  pfa -> log crp


reg <- lmer( logcrp_plasma ~ logomega3_noala_pfa + X1 + logomega3_noala_pfa*X1 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_pcoa)
summary(reg)
a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
anova(reg)[nrow(anova(reg)),ncol(anova(reg))]


ggplot(meta_pcoa, aes(X1, logomega3_noala_pfa)) + 
  geom_point(aes(color = logcrp_plasma),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff", "#b8de29ff","#29af7fff", "#2d708eff", "#440154ff", "#440154ff"))+
  labs( title = "color by logcrp_plasma",
        subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
  xlab("PCo 1")+
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


a<-quantile(meta_pcoa$logcrp_plasma, probs = seq(0, 1, 0.2))
a[1] <- -Inf
a[6] <- +Inf
meta_pcoa$logcrp_plasmaq <- cut(meta_pcoa$logcrp_plasma, breaks=a, label=c(1, 2, 3, 4, 5))

ggplot(meta_pcoa, aes(x = logcrp_plasmaq, y = logomega3_noala_pfa, fill=logcrp_plasmaq)) +
  geom_boxplot(colour = "black",lwd=0.5, outlier.color = "black",outlier.size = 1)+
  theme_bw() +
  # scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8", "#253494"))+
  ylab("omega3_noala_ddr")+
  xlab("logcrp_plasmaq")+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=1),
        plot.title=element_blank(),
        axis.title.y=element_text(color = "Black",size=10),
        axis.title.x=element_text(color = "Black",size=10),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10,angle = 300,hjust = 0,vjust=1))




# omega3  pfa -> log crp

reg <- lmer( logcrp_plasma ~ omega3_noala_pfa + X1 + omega3_noala_pfa*X1 
             +agemlvs+abx_ddr+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_pcoa)
summary(reg)
a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
anova(reg)[nrow(anova(reg)),ncol(anova(reg))]


ggplot(meta_pcoa, aes(X1, omega3_noala_pfa)) + 
  geom_point(aes(color = logcrp_plasma),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff", "#b8de29ff","#29af7fff", "#2d708eff", "#440154ff", "#440154ff"))+
  labs( title = "color by logcrp_plasma",
        subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
  xlab("PCo 1")+
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


a<-quantile(meta_pcoa$logcrp_plasma, probs = seq(0, 1, 0.2))
a[1] <- -Inf
a[6] <- +Inf
meta_pcoa$logcrp_plasmaq <- cut(meta_pcoa$logcrp_plasma, breaks=a, label=c(1, 2, 3, 4, 5))

ggplot(meta_pcoa, aes(x = logcrp_plasmaq, y = omega3_noala_pfa, fill=logcrp_plasmaq)) +
  geom_boxplot(colour = "black",lwd=0.5, outlier.color = "black",outlier.size = 1)+
  theme_bw() +
  # scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8", "#253494"))+
  ylab("omega3_noala_ddr")+
  xlab("logcrp_plasmaq")+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=1),
        plot.title=element_blank(),
        axis.title.y=element_text(color = "Black",size=10),
        axis.title.x=element_text(color = "Black",size=10),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10,angle = 300,hjust = 0,vjust=1))




# omega3  pfa -> log crp

reg <- lmer( logepa_pfa ~ logepa_ddr + X1 + logepa_ddr*X1 
             +agemlvs#+abx_ddr#+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_pcoa)
summary(reg)
a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
anova(reg)[nrow(anova(reg)),ncol(anova(reg))]


ggplot(meta_pcoa, aes(X1, logepa_ddr)) + 
  geom_point(aes(color = logepa_pfa),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff", "#b8de29ff","#29af7fff", "#2d708eff", "#440154ff", "#440154ff"))+
  labs( title = "color by logepa_pfa",
        subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
  xlab("PCo 1")+
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


a<-quantile(meta_pcoa$logepa_pfa, probs = seq(0, 1, 0.2))
a[1] <- -Inf
a[6] <- +Inf
meta_pcoa$logepa_pfaq <- cut(meta_pcoa$logepa_pfa, breaks=a, label=c(1, 2, 3, 4, 5))

ggplot(meta_pcoa, aes(x = logepa_pfaq, y = logepa_ddr, fill=logepa_pfaq)) +
  geom_boxplot(colour = "black",lwd=0.5, outlier.color = "black",outlier.size = 1)+
  theme_bw() +
  # scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8", "#253494"))+
  ylab("logepa_ddr")+
  xlab("logepa_pfaq")+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=1),
        plot.title=element_blank(),
        axis.title.y=element_text(color = "Black",size=10),
        axis.title.x=element_text(color = "Black",size=10),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10,angle = 300,hjust = 0,vjust=1))

#dpa

reg <- lmer( logdpa_pfa ~ logdpa_ddr + X1 + logdpa_ddr*X1 
             +agemlvs#+abx_ddr#+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_pcoa)
summary(reg)
a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
anova(reg)[nrow(anova(reg)),ncol(anova(reg))]


ggplot(meta_pcoa, aes(X1, logdpa_ddr)) + 
  geom_point(aes(color = logdpa_pfa),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff", "#b8de29ff","#29af7fff", "#2d708eff", "#440154ff", "#440154ff"))+
  labs( title = "color by logdpa_pfa",
        subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
  xlab("PCo 1")+
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


a<-quantile(meta_pcoa$logdpa_pfa, probs = seq(0, 1, 0.2))
a[1] <- -Inf
a[6] <- +Inf
meta_pcoa$logdpa_pfaq <- cut(meta_pcoa$logdpa_pfa, breaks=a, label=c(1, 2, 3, 4, 5))

ggplot(meta_pcoa, aes(x = logdpa_pfaq, y = logdpa_ddr, fill=logdpa_pfaq)) +
  geom_boxplot(colour = "black",lwd=0.5, outlier.color = "black",outlier.size = 1)+
  theme_bw() +
  # scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8", "#253494"))+
  ylab("logdpa_ddr")+
  xlab("logdpa_pfaq")+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=1),
        plot.title=element_blank(),
        axis.title.y=element_text(color = "Black",size=10),
        axis.title.x=element_text(color = "Black",size=10),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10,angle = 300,hjust = 0,vjust=1))


#dha

reg <- lmer( logdha_pfa ~ logdha_ddr + X1 + logdha_ddr*X1 
             +agemlvs#+abx_ddr#+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_pcoa)
summary(reg)
a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
anova(reg)[nrow(anova(reg)),ncol(anova(reg))]


ggplot(meta_pcoa, aes(X1, logdha_ddr)) + 
  geom_point(aes(color = logdha_pfa),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff", "#b8de29ff","#29af7fff", "#2d708eff", "#440154ff", "#440154ff"))+
  labs( title = "color by logdha_pfa",
        subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
  xlab("PCo 1")+
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


a<-quantile(meta_pcoa$logdha_pfa, probs = seq(0, 1, 0.2))
a[1] <- -Inf
a[6] <- +Inf
meta_pcoa$logdha_pfaq <- cut(meta_pcoa$logdha_pfa, breaks=a, label=c(1, 2, 3, 4, 5))

ggplot(meta_pcoa, aes(x = logdha_pfaq, y = logdha_ddr, fill=logdha_pfaq)) +
  geom_boxplot(colour = "black",lwd=0.5, outlier.color = "black",outlier.size = 1)+
  theme_bw() +
  # scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8", "#253494"))+
  ylab("logdha_ddr")+
  xlab("logdha_pfaq")+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=1),
        plot.title=element_blank(),
        axis.title.y=element_text(color = "Black",size=10),
        axis.title.x=element_text(color = "Black",size=10),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10,angle = 300,hjust = 0,vjust=1))




#ala

reg <- lmer( logala_pfa ~ logala_ddr + X1 + logala_ddr*X1 
             +agemlvs#+abx_ddr#+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_pcoa)
summary(reg)
a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
anova(reg)[nrow(anova(reg)),ncol(anova(reg))]


ggplot(meta_pcoa, aes(X1, logala_ddr)) + 
  geom_point(aes(color = logala_pfa),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff", "#b8de29ff","#29af7fff", "#2d708eff", "#440154ff", "#440154ff"))+
  labs( title = "color by logala_pfa",
        subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
  xlab("PCo 1")+
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


a<-quantile(meta_pcoa$logala_pfa, probs = seq(0, 1, 0.2))
a[1] <- -Inf
a[6] <- +Inf
meta_pcoa$logala_pfaq <- cut(meta_pcoa$logala_pfa, breaks=a, label=c(1, 2, 3, 4, 5))

ggplot(meta_pcoa, aes(x = logala_pfaq, y = logala_ddr, fill=logala_pfaq)) +
  geom_boxplot(colour = "black",lwd=0.5, outlier.color = "black",outlier.size = 1)+
  theme_bw() +
  # scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8", "#253494"))+
  ylab("logala_ddr")+
  xlab("logala_pfaq")+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=1),
        plot.title=element_blank(),
        axis.title.y=element_text(color = "Black",size=10),
        axis.title.x=element_text(color = "Black",size=10),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10,angle = 300,hjust = 0,vjust=1))



#tot omega3

reg <- lmer( logomega3_pfa ~ logomega3_ddr + X1 + logomega3_ddr*X1 
             #+agemlvs#+abx_ddr#+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_pcoa)
summary(reg)
a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
anova(reg)[nrow(anova(reg)),ncol(anova(reg))]


ggplot(meta_pcoa, aes(X1, logomega3_ddr)) + 
  geom_point(aes(color = logomega3_pfa),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff", "#b8de29ff","#29af7fff", "#2d708eff", "#440154ff", "#440154ff"))+
  labs( title = "color by logomega3_pfa",
        subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
  xlab("PCo 1")+
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


a<-quantile(meta_pcoa$logomega3_pfa, probs = seq(0, 1, 0.2))
a[1] <- -Inf
a[6] <- +Inf
meta_pcoa$logomega3_pfaq <- cut(meta_pcoa$logomega3_pfa, breaks=a, label=c(1, 2, 3, 4, 5))

ggplot(meta_pcoa, aes(x = logomega3_pfaq, y = logomega3_ddr, fill=logomega3_pfaq)) +
  geom_boxplot(colour = "black",lwd=0.5, outlier.color = "black",outlier.size = 1)+
  theme_bw() +
  # scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8", "#253494"))+
  ylab("logomega3_ddr")+
  xlab("logomega3_pfaq")+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=1),
        plot.title=element_blank(),
        axis.title.y=element_text(color = "Black",size=10),
        axis.title.x=element_text(color = "Black",size=10),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10,angle = 300,hjust = 0,vjust=1))



# log AA

reg <- lmer( logAA_pfa ~ logomega3_ddr + X1 + logomega3_ddr*X1 
             #+agemlvs#+abx_ddr#+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_pcoa)
summary(reg)
a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
anova(reg)[nrow(anova(reg)),ncol(anova(reg))]


ggplot(meta_pcoa, aes(X1, logomega3_ddr)) + 
  geom_point(aes(color = logAA_pfa),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff", "#b8de29ff","#29af7fff", "#2d708eff", "#440154ff", "#440154ff"))+
  labs( title = "color by logAA_pfa",
        subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
  xlab("PCo 1")+
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


a<-quantile(meta_pcoa$logAA_pfa, probs = seq(0, 1, 0.2))
a[1] <- -Inf
a[6] <- +Inf
meta_pcoa$logAA_pfaq <- cut(meta_pcoa$logAA_pfa, breaks=a, label=c(1, 2, 3, 4, 5))

ggplot(meta_pcoa, aes(x = logAA_pfaq, y = logomega3_ddr, fill=logAA_pfaq)) +
  geom_boxplot(colour = "black",lwd=0.5, outlier.color = "black",outlier.size = 1)+
  theme_bw() +
  # scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8", "#253494"))+
  ylab("logomega3_ddr")+
  xlab("logAA_pfaq")+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=1),
        plot.title=element_blank(),
        axis.title.y=element_text(color = "Black",size=10),
        axis.title.x=element_text(color = "Black",size=10),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10,angle = 300,hjust = 0,vjust=1))

#### AA omega3 ####
reg <- lmer( AA_pfa ~ logomega3_ddr + X1 + logomega3_ddr*X1 
             #+agemlvs#+abx_ddr#+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_pcoa)
summary(reg)
a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
anova(reg)[nrow(anova(reg)),ncol(anova(reg))]


ggplot(meta_pcoa, aes(X1, logomega3_ddr)) + 
  geom_point(aes(color = AA_pfa),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff", "#b8de29ff","#29af7fff", "#2d708eff", "#440154ff", "#440154ff"))+
  labs( title = "color by AA_pfa",
        subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
  xlab("PCo 1")+
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


a<-quantile(meta_pcoa$AA_pfa, probs = seq(0, 1, 0.2))
a[1] <- -Inf
a[6] <- +Inf
meta_pcoa$AA_pfaq <- cut(meta_pcoa$AA_pfa, breaks=a, label=c(1, 2, 3, 4, 5))

ggplot(meta_pcoa, aes(x = AA_pfaq, y = logomega3_ddr, fill=AA_pfaq)) +
  geom_boxplot(colour = "black",lwd=0.5, outlier.color = "black",outlier.size = 1)+
  theme_bw() +
  # scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8", "#253494"))+
  ylab("logomega3_ddr")+
  xlab("AA_pfaq")+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=1),
        plot.title=element_blank(),
        axis.title.y=element_text(color = "Black",size=10),
        axis.title.x=element_text(color = "Black",size=10),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10,angle = 300,hjust = 0,vjust=1))

#### AA omega3 noala ####
reg <- lmer( AA_pfa ~ logomega3_noala_ddr + X1 + logomega3_noala_ddr*X1 
             #+agemlvs#+abx_ddr#+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_pcoa)
summary(reg)
a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
anova(reg)[nrow(anova(reg)),ncol(anova(reg))]


ggplot(meta_pcoa, aes(X1, logomega3_noala_ddr)) + 
  geom_point(aes(color = AA_pfa),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff", "#b8de29ff","#29af7fff", "#2d708eff", "#440154ff", "#440154ff"))+
  labs( title = "color by AA_pfa",
        subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
  xlab("PCo 1")+
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


a<-quantile(meta_pcoa$AA_pfa, probs = seq(0, 1, 0.2))
a[1] <- -Inf
a[6] <- +Inf
meta_pcoa$AA_pfaq <- cut(meta_pcoa$AA_pfa, breaks=a, label=c(1, 2, 3, 4, 5))

ggplot(meta_pcoa, aes(x = AA_pfaq, y = logomega3_noala_ddr, fill=AA_pfaq)) +
  geom_boxplot(colour = "black",lwd=0.5, outlier.color = "black",outlier.size = 1)+
  theme_bw() +
  # scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8", "#253494"))+
  ylab("logomega3_noala_ddr")+
  xlab("AA_pfaq")+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=1),
        plot.title=element_blank(),
        axis.title.y=element_text(color = "Black",size=10),
        axis.title.x=element_text(color = "Black",size=10),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10,angle = 300,hjust = 0,vjust=1))

#### AA omega3 noala ####
reg <- lmer( AA_pfa ~ logomega310v + X1 + logomega310v*X1 
             #+agemlvs#+abx_ddr#+probx_ddr+bristol_ddr+ smoke12 +calor_fs_dr_ddr 
             +(1 | ID1), data=meta_pcoa)
summary(reg)
a<- cbind(beta = reg@beta[length(reg@beta)],  anova(reg)[nrow(anova(reg)),6])
anova(reg)[nrow(anova(reg)),ncol(anova(reg))]


ggplot(meta_pcoa, aes(X1, logomega310v)) + 
  geom_point(aes(color = AA_pfa),size=2, alpha=0.8, position = position_jitter()) + 
  scale_color_gradientn(colours = c("#fde725ff", "#b8de29ff","#29af7fff", "#2d708eff", "#440154ff", "#440154ff"))+
  labs( title = "color by AA_pfa",
        subtitle = paste( "P_int = ", round (anova(reg)[nrow(anova(reg)),ncol(anova(reg))], 5)))+
  xlab("PCo 1")+
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


a<-quantile(meta_pcoa$AA_pfa, probs = seq(0, 1, 0.2))
a[1] <- -Inf
a[6] <- +Inf
meta_pcoa$AA_pfaq <- cut(meta_pcoa$AA_pfa, breaks=a, label=c(1, 2, 3, 4, 5))

ggplot(meta_pcoa, aes(x = AA_pfaq, y = logomega310v, fill=AA_pfaq)) +
  geom_boxplot(colour = "black",lwd=0.5, outlier.color = "black",outlier.size = 1)+
  theme_bw() +
  # scale_fill_gradientn(colours = c("#ffffcc","#a1dab4","#41b6c4","#2c7fb8", "#253494"))+
  ylab("logomega310v")+
  xlab("AA_pfaq")+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black",size=1),
        plot.title=element_blank(),
        axis.title.y=element_text(color = "Black",size=10),
        axis.title.x=element_text(color = "Black",size=10),
        axis.text.y=element_text(color = "Black",size=10),
        axis.text.x=element_text(color = "Black",size=10,angle = 300,hjust = 0,vjust=1))



