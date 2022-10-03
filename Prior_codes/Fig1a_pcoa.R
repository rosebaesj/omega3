
rm(list=ls())
getwd()
setwd("/noGit")
options(max.print=1000000)

library(tidyverse)
library(vegan)
library(RColorBrewer)
library(grid)
library(ggpubr)
library(viridis)

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


# read in metadata

meta_stn <- data.frame(read.table("meta_stn.tsv", header=TRUE, sep='\t', check.names=F))
meta_stn$omega3_noala10v

colnames(meta_stn)
med_score<-subset(meta_stn, select = c(id_st, omega3_noala10v)) ## omega3 noala
pcoap<-inner_join(pcoap, med_score, by="id_st")

species_filt$id_st<-rownames(species_filt)
ttax_pcoap<-inner_join(pcoap, species_filt, by="id_st")
write.csv(ttax_pcoap, file="./pcoa_species2.csv")
write.csv(pcoap, file="./pcoa_species.csv")

# PCOA figure colored by act_paqq
gg<- 
  ggplot(data = pcoap, aes(X1, X2, colour = log(pcoap$omega3_noala10v)
)) +
  geom_point(size=2, alpha=0.7) + 
  xlab(paste('PCo1 (', round(100*mod$eig[1]/sum(mod$eig), 2), '%)', sep="")) + 
    ylab(paste('PCo2 (', round(100*mod$eig[2]/sum(mod$eig), 2), '%)', sep="")) +
  theme_classic()+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.line = element_line(color = "Black", size=0.5)
  ) +
  scale_color_viridis(name = "omega3 intake (no ALA)")



png(file="./figure_generated/pcoa_omega2_noala.png",
    width=5,height=5, units="in", res = 1000)
grid.arrange(gg, ncol=1, nrow=1)
dev.off()



meta_species<-read.csv(file="./data_generated/meta_species.csv", header = TRUE)
summary(meta_species$paee_pamwk)

# PCOA figure colored by physical activity: act_paqlong
#"#360167","#FF8E15",  "#FFCF07","#57c84d","#3c8321","#307e22", "#257a24""",
#"#98ff98","#7be27d","#5fc663","#41aa4a","#005b00","#004200","#002d00","#000d00"
png(file = "./figure_generated/pcoa_by_act_paqlong.png",width = 2400,height = 2400,pointsize = 50)
ggplot(data = meta_species, aes(X1, X2)) + geom_point(aes(color = act_paqlong),size=22, alpha=0.8, position = position_jitter()) +
  scale_color_gradientn(colours = c("#FDE725FF", "#B4DD2CFF","#5DC963FF","#21908CFF","#3B528BFF","#482878FF","#440154FF"))+
  ylab('PCo2 (8.8%)') + xlab('PCo1 (9.5%)') +
  theme_classic()+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 80),
        legend.text = element_text(size = 70),
        axis.title.x = element_text(size = 80),
        axis.title.y = element_text(size = 80),
        axis.text.y = element_text(size = 80),
        axis.text.x = element_text(size = 80),
        axis.line = element_line(color = "Black", size=3)
  )
dev.off()

# PCOA figure colored by physical activity: paee_pam
png(file = "./figure_generated/pcoa_by_paee_pamwk.png",width = 2400,height = 2400,pointsize = 50)
ggplot(data = meta_species, aes(X1, X2)) + geom_point(aes(color = paee_pamwk),size=22, alpha=0.8, position = position_jitter()) +
  scale_color_gradientn(colours = c("#FDE725FF", "#B4DD2CFF","#5DC963FF","#21908CFF","#3B528BFF","#482878FF","#440154FF"))+
  ylab('PCo2 (8.8%)') + xlab('PCo1 (9.5%)') +
  theme_classic()+
  theme(legend.position = "bottom",
        legend.title = element_text(size = 80),
        legend.text = element_text(size = 70),
        axis.title.x = element_text(size = 80),
        axis.title.y = element_text(size = 80),
        axis.text.y = element_text(size = 80),
        axis.text.x = element_text(size = 80),
        axis.line = element_line(color = "Black", size=3)
  )
dev.off()





## try the sam thing with not filtered data
all_species2 <- t(all_species1/100)
  
# calculate bray-crutis matrix
all_species2_bc <- vegdist(all_species2, "bray")

# PCOA using cmdscale
mm <- cmdscale(all_species2_bc, eig = T)
pcoap <- data.frame(mm$points)
pcoap$id_st<-rownames(pcoap)

# percentages of variation explained by PCO1 & 2
mm$eig[1]/sum(mm$eig)
mm$eig[2]/sum(mm$eig)

####***still not the same plz check*****##
