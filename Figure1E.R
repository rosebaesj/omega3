# PCoA

library(phyloseq)
library(ggpmisc)
library(entropart)
library(ggplot2)
library(stringr)
library(vegan)

theme_set(
  theme_bw()+
    theme(axis.line = element_line(size=1),
          axis.ticks = element_line(size=1),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face="bold"),
          plot.subtitle = element_text(hjust = 0.5)
    ))

theme_set(
  theme_bw()+
    nature_theme
)

#not do this just use species_all

species_all #this is type of OTU table


tOTU <- data.frame(t(OTU))

tOTU$id <- substr(rownames(tOTU), start = 1, stop = 6)


tOTU <- sapply(tOTU[,1:ncol(tOTU)], as.numeric)
tOTU <- as.data.frame(tOTU)
tOTUmean <- aggregate(.~id, tOTU, mean)
# gutmean <- rename(gutmean,ID=id)

rownames(tOTUmean) = tOTUmean$id


otu <- otu_table(t(tOTUmean), taxa_are_rows = TRUE)


##metadata

META <- data.frame(z_mlvs_exposure_plasma)
rownames(META) <- META$id
meta <- sample_data(META)

## tax

TAX
tax <- tax_table(TAX)





#Make Phyloseq Object



taxa_names(tax)
taxa_names(otu)

sample_names(otu)
sample_names(meta)


phy <- phyloseq(otu, tax, meta)




# #Filter
# library(metagMisc)
# phy <- phyloseq_filter_prevalence(phy, prev.trh = 0.10, abund.trh = 0.01, threshold_condition = "AND") 
# gut = t(phy@otu_table@.Data)
# dim(gut)

#Alpha Diversity
Shannon <- estimate_richness(phy, measures = "Shannon")
Simpson <- estimate_richness(phy, measures = "Simpson")
InvSimpson <- estimate_richness(phy, measures = "InvSimpson")

#Beta Diversity
Shannon_beta <- diversity(tOTUmean, index="shannon")


#richness <- data.frame(estimate_richness(phy)) ##doesn't work
###warning appears and only shannon and inv simpson were measured
Diversity <- merge(Shannon, META, by = 0, all.x = TRUE)
rownames(Diversity) <-Diversity[,1]
Diversity <- Diversity[,-1]

Diversity <- merge(Simpson, Diversity, by = 0, all.x = TRUE)
rownames(Diversity) <-Diversity[,1]
Diversity <- Diversity[,-1]

Diversity <- merge(InvSimpson, Diversity, by = 0, all.x = TRUE)
rownames(Diversity) <-Diversity[,1]
Diversity <- Diversity[,-1]


Diversity <- merge(Shannon_beta, Diversity, by = 0, all.x = TRUE)
rownames(Diversity) <-Diversity[,1]
Diversity <- Diversity[,-1]
colnames(Diversity)[1] <- c("Shannon_beta")





cor.test(Shannon$Shannon, Shannon$ala_pfa,
         alternative = c("two.sided"),
         method = c("spearman"), exact=F)
# 
# regression <- function(data = data, x = x, y = y){
#   
#   
#   x<-substitute(x)
#   y<-substitute(y)
#   l <- summary(lm(data$substitute(x)~data$substitute(x), data = data))
#   f <- str_c("noGit/result/figures/alpha_div/", deparse(x), "_", deparse(y), ".png")
#   g <- ggplot(data = data, aes(x=x, y=y))+
#     geom_point()+
#     stat_smooth(method=glm, show.legend = TRUE)+
#     labs(caption = str_c("R square = ", l$r.squared))
#   ggsave(f, plot = g, width=6, height=5, units="in", device = "png")
# 
#   }
# 
# 
# substitute(ala_pfa)
# regression(data = richness_M, x = ala_pfa, y = Shannon)


options("scipen" = 100)


##### Beta diversity Shannon ######
ggplot(data = Diversity, aes(x = , y=Shannon_beta))+
  geom_boxplot(color = "navy", outlier.alpha = 0)+
  geom_jitter(alpha = 0.3, width = 0.3,color = "navy")+
  nature_theme +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
        )


colnames(Diversity)




##### Alpha diversity ####
ggplot(data = Diversity, aes(x = "Shannon",y=Shannon))+
  geom_boxplot(color = "navy", outlier.alpha = 0)+
  geom_jitter(alpha = 0.3, width = 0.3,color = "navy")+
  nature_theme +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank())

ggsave("noGit/result/figures/alpha_div/Sahnnon_bar.png", width=3, height=2.5, units="in", device = "png")


ggplot(data = Diversity, aes(x = "Simpson",y=Simpson))+
  geom_boxplot(color = "navy", outlier.alpha = 0)+
  geom_jitter(alpha = 0.3, width = 0.3,color = "navy")+
  nature_theme +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave("noGit/result/figures/alpha_div/Simpson_bar.png", width=3, height=2.5, units="in", device = "png")


ggplot(data = Diversity, aes(x = "InvSimpson",y=InvSimpson))+
  geom_boxplot(color = "navy", outlier.alpha = 0)+
  geom_jitter(alpha = 0.3, width = 0.3,color = "navy")+
  nature_theme +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggsave("noGit/result/figures/alpha_div/InvSimpson_bar.png", width=3, height=2.5, units="in", device = "png")


# Short term intake 

ggplot(data = Diversity, aes(x=a_omega3_noala_fs_dr_w1avg, y=Shannon))+
  geom_point(alpha = 0.3, width = 0.3,color = "navy")+
  stat_smooth(method=glm)+
  labs(caption = str_c("R square = ", 
                       summary(lm(a_omega3_noala_fs_dr_w1avg~Shannon, data = Diversity))$r.squared,
                       ",\n P value = ",
                       summary(lm(crp_avg ~Shannon, data = Diversity))$coefficients[2,4]))
ggsave("noGit/result/figures/alpha_div/a_omega3_noala_fs_dr_w1avg_Sahnnon.png", width=3, height=3, units="in", device = "png")


ggplot(data = Diversity, aes(x=a_omega3_noala_fs_dr_w1avg, y=Simpson))+
  geom_point(alpha = 0.3, width = 0.3,color = "navy")+
  stat_smooth(method=glm)+
  labs(caption = str_c("R square = ", 
                       summary(lm(a_omega3_noala_fs_dr_w1avg~Simpson, data = Diversity))$r.squared,
                       ",\n P value = ",
                       summary(lm(crp_avg ~Simpson, data = Diversity))$coefficients[2,4]))
ggsave("noGit/result/figures/alpha_div/a_omega3_noala_fs_dr_w1avg_Simpson.png", width=3, height=3, units="in", device = "png")


ggplot(data = Diversity, aes(x=a_omega3_noala_fs_dr_w1avg, y=InvSimpson))+
  geom_point(alpha = 0.3, width = 0.3,color = "navy")+
  stat_smooth(method=glm)+
  labs(caption = str_c("R square = ", 
                       summary(lm(a_omega3_noala_fs_dr_w1avg~InvSimpson, data = Diversity))$r.squared,
                       ",\n P value = ",
                       summary(lm(crp_avg ~InvSimpson, data = Diversity))$coefficients[2,4]))
ggsave("noGit/result/figures/alpha_div/a_omega3_noala_fs_dr_w1avg_InvSimpson.png", width=3, height=3, units="in", device = "png")

## all of them no significant 





ggplot(data = richness_M, aes(x=ala_pfa, y=Shannon))+
  geom_point()+
  stat_smooth(method=glm)+
  labs(caption = str_c("R square = ", 
                       summary(lm(ala_pfa~Shannon, data = richness_M))$r.squared))
ggsave("noGit/result/figures/alpha_div/ala_pfa_Sahnnon.png", width=6, height=5, units="in", device = "png")

ggplot(data = richness_M, aes(x=epa_pfa, y=Shannon))+
  geom_point()+
  stat_smooth(method=glm)+
  labs(caption = str_c("R square = ", 
                       summary(lm(epa_pfa~Shannon, data = richness_M))$r.squared))
ggsave("noGit/result/figures/alpha_div/epa_pfa_Sahnnon.png", width=6, height=5, units="in", device = "png")

ggplot(data = richness_M, aes(x=dha_pfa, y=Shannon))+
  geom_point()+
  stat_smooth(method=glm)+
  labs(caption = str_c("R square = ", 
                       summary(lm(dha_pfa~Shannon, data = richness_M))$r.squared))
ggsave("noGit/result/figures/alpha_div/dha_pfa_Sahnnon.png", width=6, height=5, units="in", device = "png")

ggplot(data = richness_M, aes(x=dpa_pfa, y=Shannon))+
  geom_point()+
  stat_smooth(method=glm)+
  labs(caption = str_c("R square = ", 
                       summary(lm(dpa_pfa~Shannon, data = richness_M))$r.squared))
ggsave("noGit/result/figures/alpha_div/dpa_pfa_Sahnnon.png", width=6, height=5, units="in", device = "png")






ggplot(data = richness_M, aes(x=a_omega3_fs_dr_w1avg, y=Shannon))+
  geom_point()+
  stat_smooth(method=glm)+
  labs(caption = str_c("R square = ", 
                       summary(lm(a_omega3_fs_dr_w1avg~Shannon, data = richness_M))$r.squared))
ggsave("noGit/result/figures/alpha_div/a_omega3_fs_dr_w1avg_Sahnnon.png", width=6, height=5, units="in", device = "png")

ggplot(data = richness_M, aes(x=a_omega3_noala_fs_dr_w1avg, y=Shannon))+
  geom_point()+
  stat_smooth(method=glm)+
  labs(caption = str_c("R square = ", 
                       summary(lm(a_omega3_noala_fs_dr_w1avg~Shannon, data = richness_M))$r.squared))
ggsave("noGit/result/figures/alpha_div/a_omega3_noala_fs_dr_w1avg_Sahnnon.png", width=6, height=5, units="in", device = "png")

ggplot(data = richness_M, aes(x=omega3_ffq1, y=Shannon))+
  geom_point()+
  stat_smooth(method=glm)+
  labs(caption = str_c("R square = ", 
                       summary(lm(omega3_ffq1~Shannon, data = richness_M))$r.squared))
ggsave("noGit/result/figures/alpha_div/omega3_ffq1_Sahnnon.png", width=3, height=3, units="in", device = "png")


ggplot(data = richness_M, aes(x=omega3_noala_ffq1, y=Shannon))+
  geom_point()+
  stat_smooth(method=glm)+
  labs(caption = str_c("R square = ", 
                       summary(lm(omega3_noala_ffq1 ~Shannon, data = richness_M))$r.squared))
ggsave("noGit/result/figures/alpha_div/omega3_noala_ffq1_Sahnnon.png", width=6, height=5, units="in", device = "png")

ggplot(data = richness_M, aes(x=crp_avg, y=Shannon))+
  geom_point()+
  stat_smooth(method=glm)+
  labs(caption = str_c("R square = ", 
                       summary(lm(crp_avg ~Shannon, data = richness_M))$r.squared, 
                       ", P value = ",
                      summary(lm(crp_avg ~Shannon, data = richness_M))$coefficients[2,4]))
ggsave("noGit/result/figures/alpha_div/crp_avg_Sahnnon.png", width=6, height=5, units="in", device = "png")






##### intakes and plasma #####

ggplot(data = richness_M, aes(x=epa_ffq1, y=epa_pfa))+
  geom_point()+
  stat_smooth(method=glm)+
  labs(caption = str_c("R square = ", 
                       summary(lm(epa_ffq1 ~epa_pfa, data = richness_M))$r.squared,
                       ",\n P value = ",
                       summary(lm(epa_ffq1 ~epa_pfa, data = richness_M))$coefficients[2,4]))
ggsave("noGit/result/figures/metabolites/epa_ffa1_pfa.png", width=3, height=3, units="in", device = "png")

ggplot(data = richness_M, aes(x=ala_ffq1, y=ala_pfa))+
  geom_point()+
  stat_smooth(method=glm)+
  labs(caption = str_c("R square = ", 
                       summary(lm(ala_ffq1 ~ala_pfa, data = richness_M))$r.squared,
                       ",\n P value = ",
                       summary(lm(ala_ffq1 ~ala_pfa, data = richness_M))$coefficients[2,4]))
ggsave("noGit/result/figures/metabolites/ala_ffa1_pfa.png", width=3, height=3, units="in", device = "png")


ggplot(data = richness_M, aes(x=dpa_ffq1, y=dpa_pfa))+
  geom_point()+
  stat_smooth(method=glm)+
  labs(caption = str_c("R square = ", 
                       summary(lm(dpa_ffq1 ~dpa_pfa, data = richness_M))$r.squared,
                       ",\n P value = ",
                       summary(lm(dpa_ffq1 ~dpa_pfa, data = richness_M))$coefficients[2,4]))
ggsave("noGit/result/figures/metabolites/dpa_ffa1_pfa.png", width=3, height=3,units="in", device = "png")

ggplot(data = richness_M, aes(x=dha_ffq1, y=dha_pfa))+
  geom_point()+
  stat_smooth(method=glm)+
  labs(caption = str_c("R square = ", 
                       summary(lm(dha_ffq1 ~dha_pfa, data = richness_M))$r.squared,
                       ",\n P value = ",
                       summary(lm(dha_ffq1 ~dha_pfa, data = richness_M))$coefficients[2,4]))
ggsave("noGit/result/figures/metabolites/dha_ffa1_pfa.png", width=3, height=3, units="in", device = "png")


cor <- richness





#### Beta Diversity #####

phy_bray <- phyloseq::distance(phy, method = "bray")
ordination = ordinate(phy, method="PCoA", distance=phy_bray)
plot_ordination(phy, ordination, color="SampleType") + theme(aspect.ratio=1)


#Beta Diversity



omega3 <- data.frame(c(
                       'trans_avg', 'omega3_avg', 'omega6_avg', 'ala_avg', 'epa_avg', 'dpa_avg', 'dha_avg', 'omega3_noala_avg',  
                    'a_omega3_fs_dr_w1avg', 'a_omega3_noala_fs_dr_w1avg', 'a_trn07_fo_dr_w1avg',
                    'a_ala_fs_dr_w1avg', 'a_f205_fs_dr_w1avg', 'a_f226_fs_dr_w1avg', 'a_p22_5_fs_dr_w1avg',
                    'ala_pfa', 'epa_pfa', 'dpa_pfa', 'dha_pfa',
                    'aofib_avg', 'crp_avg', 'abx_avg', 'agemlvs', 'calor10n', 'act10', 'bmi10', 'bristol_avg'
                    ))
                    
permanova = data.frame()

for (i in 1: nrow(omega3)){
  a<- omega3[i,]
  p <- adonis2(phy_bray ~ get(a), data = META, permutations=999, method = "bray", na.action	= na.exclude)
  permanova <- rbind( permanova, new = p[1,])
  rownames(permanova)[i]<-c(a)
}

write.table

  
  
  
  
  
#Run Maaslin2
library(Maaslin2)
fit_data <- Maaslin2(
  gut, dat, 'demo_output', transform = "none", heatmap_first_n = 100,
  max_significance = 0.25,
  fixed_effects = c('cod_liv_oil_avg','linseed_avg','flax_avg','omega_3_epa_avg', 'walnuts_avg', 'avocado_avg', 'leafvegavg', 'fishavg', 'agemlvs','calor_avg','abx_avg'),
  random_effects = 'ID1',
  normalization = 'NONE',
  standardize = FALSE)


#FFQ Data
ffq1 <- read_sas("/Users/jorickbater/Desktop/Mingyang Stuff/Data/mlvsffq1_excl_nts.sas7bdat")
om1 <- ffq1[,c("id","pfa183n3c07_fs_ffq1", "pf205n3c07_fs_ffq1", "pf225n3c07_fs_ffq1", "pf226n3c07_fs_ffq1")]
ffq2 <- read_sas("/Users/jorickbater/Desktop/Mingyang Stuff/Data/mlvsffq2_excl_nts.sas7bdat")
om2 <- ffq2[,c("id","pfa183n3c07_fs_ffq2", "pf205n3c07_fs_ffq2", "pf225n3c07_fs_ffq2", "pf226n3c07_fs_ffq2")]
dr1 <- read_sas("/Users/jorickbater/Desktop/Mingyang Stuff/Data/dr_nts_wk1_mean_e_adjusted.sas7bdat")
dom1 <- dr1[,c("id","a_omega3_fs_dr_w1avg", "a_p22_5_fs_dr_w1avg")]
dr2 <- read_sas("/Users/jorickbater/Desktop/Mingyang Stuff/Data/dr_nts_wk2_mean_e_adjusted.sas7bdat")
dom2 <- dr2[,c("id","a_omega3_fs_dr_w2avg", "a_p22_5_fs_dr_w2avg")]

dat_list = list(om1,om2,dom1,dom2)

my_merge <- function(df1, df2){                               
  merge(df1, df2, by = "id")
}

totom <- Reduce(my_merge,dat_list)
totom$pfa183n3c07_fs_avg <- rowMeans(totom[,c('pfa183n3c07_fs_ffq1', 'pfa183n3c07_fs_ffq2')], na.rm=TRUE)
totom$pf205n3c07_fs_avg <- rowMeans(totom[,c('pf205n3c07_fs_ffq1', 'pf205n3c07_fs_ffq2')], na.rm=TRUE)
totom$pf225n3c07_fs_avg <- rowMeans(totom[,c('pf225n3c07_fs_ffq1', 'pf225n3c07_fs_ffq2')], na.rm=TRUE)
totom$pf226n3c07_fs_avg <- rowMeans(totom[,c('pf226n3c07_fs_ffq1', 'pf226n3c07_fs_ffq2')], na.rm=TRUE)











# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### extra code ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
met <- read.csv('noGit/mlvs_exposure_wenjie.csv', header=T, row.names = 1)
met$id = rownames(met)

z_metadata_all$choc_avg <- (z_metadata_all$dchocffq1+z_metadata_all$mchocffq1+z_metadata_all$dchocffq2+z_metadata_all$mchocffq2)/2
z_metadata_all$choc_avg[is.na(z_metadata_all$choc_avg)] <- 0


#Blood Fatty acids
fa <- read.table("noGit/matching_plasma_id.tsv")
fa <- as.data.frame(fa)
fa$ID <- as.numeric(fa$ID)
fa <- fa[,3:47]
dat <- merge(fa,tot,by="ID")




#below works fine but...
gut <-read_tsv("noGit/metaphlan_taxonomic_profiles.tsv")



gut <- t(gut)
cols <- gut[1,]
colnames(gut) <- cols
gut = gut[-1,]
#gut <- gut[ , grepl( "s__" , colnames(gut) ) ]
#colnames(gut) <- gsub('.*s__', '', colnames(gut))

gut <- as.data.frame(gut)
gut$id <- substr(rownames(gut), start = 1, stop = 6)


gut <- sapply(gut[,1:ncol(gut)], as.numeric)
gut <- as.data.frame(gut)
gutmean <- aggregate(.~id, gut, mean)
# gutmean <- rename(gutmean,ID=id)

rownames(gutmean) = gutmean$id





z_mlvs_exposure <- read.csv('/Users/jorickbater/Desktop/Mingyang Stuff/Data/mlvs_exposure_wenjie.csv', header=T, row.names = 1)
z_mlvs_exposure$id = rownames(z_mlvs_exposure)
z_mlvs_exposure$id <- as.factor(z_mlvs_exposure$id)
z_mlvs_exposure <- merge(totom,z_mlvs_exposure,by="id")

z_metadata_all$participant <- rownames(z_metadata_all)
met <- z_metadata_all

#Merge ID
ID <- read_sas("noGit/mlvs_bld_phase12.sas7bdat")
ID = ID[,c("HarvardID","id")]

tot <- merge(met,ID, by="id")
tot$id <- NULL
tot <- rename(tot,id=HarvardID)






#Gut Taxa
tax_rpk_name <- read_tsv("noGit/metaphlan_taxonomic_profiles.tsv")
names(tax_rpk_name)[names(tax_rpk_name) == '# taxonomy'] <- 'Sample'
tax_rpk_name<-tax_rpk_name %>%
  separate(Sample, c("kingdom",       "phylum",        "class" ,        "order",         "family",        "genus" ,        "species" ,      "strain"), 
           sep = '\\|', remove = TRUE)






# only keep species-level features
tax_rpk_species <- subset(tax_rpk_name,!is.na(species) & is.na(strain))
rownames(tax_rpk_species)<-tax_rpk_species$species
tax_rpk_species<-tax_rpk_species[,-c(1:8)]

all_species <- as.data.frame(tax_rpk_species)
all_species <- all_species[ , !colnames(all_species) %in% c(grep("15074981", colnames(all_species), value = T))]
dim(all_species)

dim(all_species[apply(all_species, 1, function(bug) sum(bug >= 0.01) >= 0.1*ncol(all_species)), ])
all_species_filt <- all_species[ apply(all_species, 1, function(bug) sum(bug >= 0.01) >= 0.1*ncol(all_species)), ]
species_list <- rownames(all_species_filt)
species_list.format <- tolower(gsub(pattern = "s__", "", species_list))

z_idkey <- read_sas("noGit/idlinks.sas7bdat")
z_idkey <- z_idkey[2:909,]
rownames(z_idkey) <- z_idkey$id
z_idkey <- z_idkey[order(z_idkey$id), ] 
z_idkey <- rename(z_idkey,ID1=aliasid8digits)



med_id<-subset(tot, select=id)
ttax_rpk <- as.data.frame(t(tax_rpk_species))
ttax_rpk$id<-rownames(ttax_rpk)
ttax_rpk$id <- gsub("([0-9]+)_.*", "\\1", ttax_rpk$id)

ttax_rpk<-inner_join(med_id, ttax_rpk, by="id")
rownames(ttax_rpk)<-ttax_rpk$id
ttax_rpk <- subset(ttax_rpk, select = -id)
View(ttax_rpk)

dat <- tot
nam = as.vector(dat$ID)
gut <- filter(gut, ID %in% nam)

gut$ID <- NULL
gut[,2:544] <- gut[,2:544]/100

rownames(dat) = dat$ID
dat$ID = as.factor(dat$ID)











#Filter
library(metagMisc)
phy <- phyloseq_filter_prevalence(phy, prev.trh = 0.10, abund.trh = 0.01, threshold_condition = "AND") 
gut = t(phy@otu_table@.Data)
dim(gut)

#Alpha Diversity
plot_richness(phy, x="ala", measures = c("Shannon", "InvSimpson")) + geom_boxplot()
plot_richness(phy, x="P10", measures = c("Shannon", "InvSimpson")) + geom_boxplot()
plot_richness(phy, x="P13", measures = c("Shannon", "InvSimpson")) + geom_boxplot()
plot_richness(phy, x="P14", measures = c("Shannon", "InvSimpson")) + geom_boxplot()
plot_richness(phy, x="crp_plasma1", measures = c("Shannon", "InvSimpson")) + geom_boxplot()
plot_richness(phy, x="omega_3_epa_ffq1", measures = c("Shannon", "InvSimpson")) + geom_boxplot()

#Beta Diversity
source("miseqR.R")
set.seed(1)
dat_scale <- phy %>%
  scale_reads(round = "round") 


dat_bray <- phyloseq::distance(dat_scale, method = "bray")
ordination = ordinate(phy, method="PCoA", distance=dat_bray)
plot_ordination(dat, ordination, color="SampleType") + theme(aspect.ratio=1)

#PERMANOVA
library(microbiome)
library(vegan)
rel <- microbiome::transform(phy, "compositional")

#Beta Diversity

p <- plot_landscape(rel, method = "NMDS", distance = "bray", col = "crp_avg", size = 3)
plot(p)

otu <- abundances(phy)
meta <- meta(phy)
permanova <- adonis(t(otu) ~ P3 + P10 + P13 + P14 + crp_avg,
data = dat, permutations=99, method = "bray")
permanova <- adonis(t(otu) ~ crp_avg,
                   data = meta, permutations=999, method = "bray")

permanova

coef <- coefficients(permanova)["crp_avg",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")

#Run Maaslin2
library(Maaslin2)
fit_data <- Maaslin2(
gut, dat, 'demo_output', transform = "none", heatmap_first_n = 100,
max_significance = 0.25,
fixed_effects = c('cod_liv_oil_avg','linseed_avg','flax_avg','omega_3_epa_avg', 'walnuts_avg', 'avocado_avg', 'leafvegavg', 'fishavg', 'agemlvs','calor_avg','abx_avg'),
random_effects = 'ID1',
normalization = 'NONE',
standardize = FALSE)


#FFQ Data
ffq1 <- read_sas("/Users/jorickbater/Desktop/Mingyang Stuff/Data/mlvsffq1_excl_nts.sas7bdat")
om1 <- ffq1[,c("id","pfa183n3c07_fs_ffq1", "pf205n3c07_fs_ffq1", "pf225n3c07_fs_ffq1", "pf226n3c07_fs_ffq1")]
ffq2 <- read_sas("/Users/jorickbater/Desktop/Mingyang Stuff/Data/mlvsffq2_excl_nts.sas7bdat")
om2 <- ffq2[,c("id","pfa183n3c07_fs_ffq2", "pf205n3c07_fs_ffq2", "pf225n3c07_fs_ffq2", "pf226n3c07_fs_ffq2")]
dr1 <- read_sas("/Users/jorickbater/Desktop/Mingyang Stuff/Data/dr_nts_wk1_mean_e_adjusted.sas7bdat")
dom1 <- dr1[,c("id","a_omega3_fs_dr_w1avg", "a_p22_5_fs_dr_w1avg")]
dr2 <- read_sas("/Users/jorickbater/Desktop/Mingyang Stuff/Data/dr_nts_wk2_mean_e_adjusted.sas7bdat")
dom2 <- dr2[,c("id","a_omega3_fs_dr_w2avg", "a_p22_5_fs_dr_w2avg")]

dat_list = list(om1,om2,dom1,dom2)

my_merge <- function(df1, df2){                               
  merge(df1, df2, by = "id")
}

totom <- Reduce(my_merge,dat_list)
totom$pfa183n3c07_fs_avg <- rowMeans(totom[,c('pfa183n3c07_fs_ffq1', 'pfa183n3c07_fs_ffq2')], na.rm=TRUE)
totom$pf205n3c07_fs_avg <- rowMeans(totom[,c('pf205n3c07_fs_ffq1', 'pf205n3c07_fs_ffq2')], na.rm=TRUE)
totom$pf225n3c07_fs_avg <- rowMeans(totom[,c('pf225n3c07_fs_ffq1', 'pf225n3c07_fs_ffq2')], na.rm=TRUE)
totom$pf226n3c07_fs_avg <- rowMeans(totom[,c('pf226n3c07_fs_ffq1', 'pf226n3c07_fs_ffq2')], na.rm=TRUE)