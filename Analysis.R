library(haven)
library(data.table)
library(dplyr)
library(tidyverse)

#setwd("/Users/jorickbater/Downloads")
getwd()
### Pre-processing steps 
source("noGit/Rstart.R")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
########################## import data ###############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# read in metadata (925)
mlvs_metadata_925 <- read.table('./noGit/maaslin2/pcl/mlvs_metadata_925.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(mlvs_metadata_925) <- mlvs_metadata_925 $sample
mlvs_metadata_925  <- mlvs_metadata_925 [,-1]
mlvs_metadata_925  <- as.data.frame(t(mlvs_metadata_925 ))

# read in species
species_all <- read.table('./noGit/maaslin2/pcl/species_all.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(species_all) <- species_all$sample
species_all <- species_all[,-1]
species_all <- as.data.frame(t(species_all))

# merge metadata and species
species_metadata <- merge(x = mlvs_metadata_925, y = species_all, by = "row.names", all = TRUE)
rownames(species_metadata) <- species_metadata[ ,1]
species_metadata[ ,1] <- NULL

ddr_plot <- species_metadata
checkit <- colnames(ddr_plot)
#View(checkit)

for(i in 1:ncol(ddr_plot)) { 
  ddr_plot[ , i] <- as.numeric(as.character(ddr_plot[, i])) 
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
####################### Fig 3A. MaAsLin2 ################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

bug_file<-c('species_all')
metadata_file <- c( 'ala_ffqcum','epa_ffqcum','dha_ffqcum', 'dpa_ffqcum', 'trans_ffqcum', 'omega6_ffqcum', 'omega3_ffqcum', 'omega3_noala_ffqcum',
                    'ala_avg_fs_ffq', 'epa_avg_fs_ffq', 'dha_avg_fs_ffq', 'dpa_avg_fs_ffq', 'trans_avg_fs_ffq', 'omega6_avg_fs_ffq', 'omega3_avg_fs_ffq', 'omega3_noala_avg_fs_ffq', 'a_omega3_fs_ddr', 'a_omega3_noala_fs_ddr', 'a_trans_fs_ddr', 'a_ala_fs_ddr', 'a_epa_fs_ddr', 'a_dha_fs_ddr', 'a_dpa_fs_ddr',
                   'biomarker_logcrp', 'biomarker_logcrp_complete', 'a_omega3_fs_ddr_logcrp', 'biomarker_loghdl', 'biomarker_loghdl_complete', 'a_omega3_fs_ddr_loghdl', 'biomarker_logtg', 'biomarker_logtg_complete', 'a_omega3_fs_ddr_logtg', 'a_trans_fs_ddr_logcrp', 'a_trans_fs_ddr_loghdl', 'a_trans_fs_ddr_logtg',
                   'a_ala_fs_ddr_logcrp', 'a_ala_fs_ddr_loghdl', 'a_ala_fs_ddr_logtg',
                   'a_epa_fs_ddr_logcrp', 'a_epa_fs_ddr_loghdl', 'a_epa_fs_ddr_logtg',
                   'a_dha_fs_ddr_logcrp', 'a_dha_fs_ddr_loghdl', 'a_dha_fs_ddr_logtg',
                   'a_dpa_fs_ddr_logcrp', 'a_dpa_fs_ddr_loghdl', 'a_dpa_fs_ddr_logtg')


# arcsin sqrt transformation  
for (a in 1:length(metadata_file)){
  # (1) with no abundance/prevalence filter since i already did this myself
  
  for (b in 1:length(bug_file))
  {
    input_data <- read.table(paste('./noGit/maaslin2/pcl/', bug_file[b], '.pcl', sep=''),
                             header = TRUE, row.names = 1
                             #  ,check.names = FALSE ##this gets rid of Xs in colnames
                             )
    input_metadata <- read.table(paste('./noGit/maaslin2/pcl/', metadata_file[a], '.pcl', sep=''),
                                 header = TRUE, row.names = 1
                                 #  ,check.names = FALSE ##this gets rid of Xs in colnames
                                 )
    Maaslin2(input_data = input_data,#= paste('./noGit/maaslin2/pcl/', bug_file[b], '.pcl', sep=''), 
             input_metadata =input_metadata, # = paste('./noGit/maaslin2/pcl/', metadata_file[a], '.pcl', sep=''),
             output           = paste('./noGit/maaslin2/omega_output_', metadata_file[a], '_', bug_file[b], '/', sep=''),
             normalization    = 'NONE', 
             standardize      = 'FALSE',
             transform        = 'AST', 
             analysis_method  = 'LM', 
             random_effects   = 'participant',
             min_abundance    = 0, 
             min_prevalence   = 0,
             cores            = 1)
    
  }
}
##takes ~ 30min

features <- ddr_plot[,c(1:36)]
microbiomes <- ddr_plot[,c(37:ncol(ddr_plot))]

Maaslin2(input_data = microbiomes, 
         input_metadata = features,
         output           = './noGit/maaslin2/microbiomes_features/',
         normalization    = 'NONE', 
         standardize      = 'FALSE',
         transform        = 'AST', 
         analysis_method  = 'LM', 
         random_effects   = 'participant',
         min_abundance    = 0, 
         min_prevalence   = 0,
         cores            = 1)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################## @@ Figure 3B/C. scatter plots ###################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

######### @ need modification from fiber to omega3# ########

# Function of positive associations
scatter_plot_positive <- function(x,y,xlab,ylab) {
  plot <- ggplot(
    data=ddr_plot, aes_string(x=x, y=y)) +
    geom_point( aes(), fill = '#D9EF8B', color = 'black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
    stat_smooth(method = "glm", color ='black')+ 
    guides(alpha='none')+labs("")+
    xlab(xlab) +  ylab(ylab)  +
    nature_theme +
    guides(legend.position=NULL)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title=element_text(size=20),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.ticks =element_blank(),
          plot.title = element_blank())
  filepath <- paste('./noGit/result/figures/scatter_species/', y, '_', x, '.png', sep ="" )
  ggsave(filename=filepath, plot=plot, width = 5, height = 5, dpi = 600) 
}# end of function

positive_bugs_for_omega3 <- c('roseburia_hominis',
                                 'bacteroides_xylanisolvens',
                                 'clostridium_sp_cag_167',
                                 'clostridium_sp_cag_253',
                                 'enterorhabdus_caecimuris',
                                 'paraprevotella_clara',
                                 'streptococcus_thermophilus',
                                 'slackia_isoflavoniconvertens',
                                 'anaerostipes_hadrus',
                                 'roseburia_intestinalis',
                                 'adlercreutzia_equolifaciens') 

for (i in 1:11){
  scatter_plot_positive('a_omega3_fs_dr',positive_bugs_for_omega3[i],"Short-term omega3\n(g/d)",
                        paste(positive_bugs_for_omega3[i], "\n(relative abundance", sep = ''))
}

scatter_plot_positive('dha10v','ruthenibacterium_lactatiformans',"dha10v\n(g/d)","ruthenibacterium_lactatiformans\n(relative abundance")
scatter_plot_positive('a_omega3_noala_fs_dr','ruthenibacterium_lactatiformans',"a_omega3_noala_fs_dr\n(g/d)","ruthenibacterium_lactatiformans\n(relative abundance")
scatter_plot_positive('a_ala_fs_dr','roseburia_hominis',"a_omega3_noala_fs_dr\n(g/d)","roseburia_hominis\n(relative abundance")
scatter_plot_positive('a_ala_fs_dr','roseburia_hominis',"a_omega3_noala_fs_dr\n(g/d)","roseburia_hominis\n(relative abundance")

scatter_plot_positive('a_trn07_fo_dr','roseburia_hominis',"a_omega3_noala_fs_dr\n(g/d)","roseburia_hominis\n(relative abundance")
scatter_plot_positive('logcrp','bacteroides_uniformis',"Log CRP","Bacteroides uniformis\n(relative abundance)")


# Function of negative associations
scatter_plot_negative <- function(x,y,xlab,ylab) {
  plot <- ggplot(
    data=ddr_plot, aes_string(x=x, y=y)) +
    geom_point( aes(), fill = '#C92D39', color = 'black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
    stat_smooth(method = "glm", color ='black')+ 
    guides(alpha='none')+labs("")+
    xlab(xlab) +  ylab(ylab)  +
    nature_theme +
    guides(legend.position=NULL)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title=element_text(size=20),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.ticks =element_blank(),
          plot.title = element_blank())
  filepath <- paste('./noGit/result/figures/scatter_species/', y, '_', x, '.png', sep ="" )
  ggsave(filename=filepath, plot=plot, width = 5, height = 5, dpi = 600) 
}# end of function

negative_bugs_for_omega3 <- c("bifidobacterium_bifidum",
                              'ruminococcaceae_bacterium_d5',
                              'dialister_invisus',
                              'clostridium_leptum',
                              'enterorhabdus_caecimuris',
                              'bifidobacterium_adolescentis',
                              'alistipes_putredinis',
                              'odoribacter_splanchnicus',
                              'collinsella_stercoris',
                              'erysipelatoclostridium_ramosum',
                              'bilophila_wadsworthia')


for (i in 1:11){
  scatter_plot_negative('a_omega3_fs_dr',negative_bugs_for_omega3[i],"Short-term omega3\n(g/d)",
                        paste(negative_bugs_for_omega3[i], "\n(relative abundance", sep = ''))
}


catter_plot_negative('logcrp','eubacterium_eligens',"Short-term dietary fiber\n(g/d)","Eubacterium eligens\n(relative abundance")



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
######### Figure 2. interaction ##***between DDR fiber and P. copri on logCRP****#####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################# read in data ###################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# read in metadata (925)
mlvs_metadata_925 <- read.table('./noGit/maaslin2/pcl/mlvs_metadata_925.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(mlvs_metadata_925) <- mlvs_metadata_925$sample
mlvs_metadata_925  <- mlvs_metadata_925 [,-1]
mlvs_metadata_925  <- as.data.frame(t(mlvs_metadata_925 ))

# read in species
species_all <- read.table('./noGit/maaslin2/pcl/species_all.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(species_all) <- species_all$sample
species_all <- species_all[,-1]
species_all <- as.data.frame(t(species_all))

# read in MDS1 and MDS2
bugs_pcoa <- species_all
bugs_pcoa <- capscale(t(bugs_pcoa) ~ 1, distance = 'bray') #Question? old pcoa used asin sqrt transformed bug data?
variance.bugs <- head(eigenvals(bugs_pcoa)/sum(eigenvals(bugs_pcoa)))
x_variance.bugs <- as.integer(variance.bugs[1]*100)
y_variance.bugs <- as.integer(variance.bugs[2]*100)

# sample scores
ord.scores <- as.data.frame(scores(bugs_pcoa, display = "sites"))

# bug scores
ord.bug.scores <- as.data.frame(scores(bugs_pcoa, display = "species") )
#View(abs(ord.bug.scores))

# merge the data (change names for trans)
association <- mlvs_metadata_925[ , colnames(mlvs_metadata_925) %in% c('participant', 'agemlvs','abx', 'bmi12', 
                                    'logcrp','loghdl','logtg','ala10v','epa10v','dha10v', 'dpa10v', 'calor10v', 
                                    'a_omega3_fs_dr', 'a_omega3_noala_fs_dr', 'calor_fs_dr','act10', 'bristolcat')] 

#merge association data with MDS bug scores
association_mds <- merge(x = association, y = ord.bug.scores, by = "row.names", all = TRUE)

# set first column as rownames and remove it
rownames(association_mds) <- association_mds[ ,1]
association_mds[ ,1] <- NULL

# merge with the bug data
mns <- colMeans(species_all, na.rm=TRUE)
order(-mns)
species_all_sort <- species_all[,order(-mns)]
colMeans(species_all_sort)

association_bugs <- merge(x = association_mds, y = species_all_sort, by = 'row.names', all = TRUE)
rownames(association_bugs) <- association_bugs[ ,1]
association_bugs[ ,1] <- NULL

head(association_bugs)
#sort by MDS1
association_bugs_sort <- association_bugs[order(association_bugs$MDS1),]

association_bugs$prevotella_copri2c <- with(association_bugs, ifelse(prevotella_copri>0, 1, 0))
table(association_bugs$prevotella_copri2c)
#0   1 
#693 232

#association_bugs$lactobacillus_rogosae <- with(association_bugs, ifelse(lactobacillus_rogosae>0, 1, 0))
#association_bugs$bifidobacterium_adolescentis <- with(association_bugs, ifelse(bifidobacterium_adolescentis>0, 1, 0))
#association_bugs$bifidobacterium_longum<- with(association_bugs, ifelse(bifidobacterium_longum>0, 1, 0))
#association_bugs$bifidobacterium_pseudocatenulatum<-- with(association_bugs, ifelse(bifidobacterium_pseudocatenulatum>0, 1, 0))
#association_bugs$bifidobacterium_bifidum<- with(association_bugs, ifelse(bifidobacterium_bifidum>0, 1, 0))

library(lmerTest)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##################### interaction model with MDS1 ####################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#crude interaction model: n=775, p-int=0.663 (1), 0.031 (2)
lm.crp <- lmer(formula = logcrp ~ a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + (1 | participant), data=association_bugs)
anova(lm.crp)
summary(lm.crp)

lm.crp <- lmer(formula = logcrp ~ a_omega3_fs_dr + MDS2 + a_omega3_fs_dr*MDS2 + (1 | participant), data=association_bugs)
anova(lm.crp)
summary(lm.crp)

#further adjusting for age: n=775, p-int=0.689 (1), 0.037 (2)
lm.crp <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + (1 | participant), data=association_bugs)
anova(lm.crp)
summary(lm.crp)

lm.crp <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + MDS2 + a_omega3_fs_dr*MDS2 + (1 | participant), data=association_bugs)
anova(lm.crp)
summary(lm.crp)

#main model, try further adjusting for some covariates: age, antibiotics, calorie; n=773, p-int=0.749 (1), 0.044 (2)
lm.crp <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + calor_fs_dr + abx + (1 | participant), data=association_bugs)
anova(lm.crp)
summary(lm.crp)

lm.crp <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + MDS2 + a_omega3_fs_dr*MDS2 + calor_fs_dr + abx + (1 | participant), data=association_bugs)
anova(lm.crp)
summary(lm.crp)

#additionally adjusting for bristolcat or physical activity did not change the results; n=752, p-int=0.9113 (1), 0.085 (2)
association_bugs$bristolcat <- as.factor(association_bugs$bristolcat)
lm.crp <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + calor_fs_dr + abx + bristolcat + act10 + (1 | participant), data=association_bugs)
anova(lm.crp)
summary(lm.crp)

lm.crp <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + MDS2 + a_omega3_fs_dr*MDS2 + calor_fs_dr + abx + bristolcat + act10 + (1 | participant), data=association_bugs)
anova(lm.crp)
summary(lm.crp)

############HDL
#crude interaction model: n=775, p-int=0.663
lm.hdl <- lmer(formula = loghdl ~ a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + (1 | participant), data=association_bugs)
anova(lm.hdl)
summary(lm.hdl)

#further adjusting for age: n=775, p-int=0.689
lm.hdl <- lmer(formula = loghdl ~ agemlvs + a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + (1 | participant), data=association_bugs)
anova(lm.hdl)
summary(lm.hdl)

#main model, try further adjusting for some covariates: age, antibiotics, calorie; n=773, p-int=0.749
lm.hdl <- lmer(formula = loghdl ~ agemlvs + a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + calor_fs_dr + abx + (1 | participant), data=association_bugs)
anova(lm.hdl)
summary(lm.hdl)

#additionally adjusting for bristolcat or physical activity did not change the results; n=752, p-int=0.9113
association_bugs$bristolcat <- as.factor(association_bugs$bristolcat)
lm.hdl <- lmer(formula = loghdl ~ agemlvs + a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + calor_fs_dr + abx + bristolcat + act10 + (1 | participant), data=association_bugs)
anova(lm.hdl)
summary(lm.hdl)

############TG############
#crude interaction model: n=775, p-int=0.663
lm.tg <- lmer(formula = logtg ~ a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + (1 | participant), data=association_bugs)
anova(lm.tg)
summary(lm.tg)

#further adjusting for age: n=775, p-int=0.689
lm.tg <- lmer(formula = logtg ~ agemlvs + a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + (1 | participant), data=association_bugs)
anova(lm.tg)
summary(lm.tg)

#main model, try further adjusting for some covariates: age, antibiotics, calorie; n=773, p-int=0.749
lm.tg <- lmer(formula = logtg ~ agemlvs + a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + calor_fs_dr + abx + (1 | participant), data=association_bugs)
anova(lm.tg)
summary(lm.tg)

#additionally adjusting for bristolcat or physical activity did not change the results; n=752, p-int=0.9113
association_bugs$bristolcat <- as.factor(association_bugs$bristolcat)
lm.tg <- lmer(formula = logtg ~ agemlvs + a_omega3_fs_dr + MDS1 + a_omega3_fs_dr*MDS1 + calor_fs_dr + abx + bristolcat + act10 + (1 | participant), data=association_bugs)
anova(lm.tg)
summary(lm.tg)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
########## interaction model with P copri ############
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#model adjusting for age only, p-int=0.49 ########
lm.crp.prevotella_copri <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + prevotella_copri2c + a_omega3_fs_dr*prevotella_copri2c + (1 | participant), data=association_bugs)
anova(lm.crp.prevotella_copri)
summary(lm.crp.prevotella_copri)

#lm.crp.lactobacillus_rogosae <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + lactobacillus_rogosae + a_omega3_fs_dr*lactobacillus_rogosae + (1 | participant), data=association_bugs)
#anova(lm.crp.lactobacillus_rogosae)
#summary(lm.crp.lactobacillus_rogosae)

#lm.crp.bifidobacterium_adolescentis <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + bifidobacterium_adolescentis + a_omega3_fs_dr*bifidobacterium_adolescentis + (1 | participant), data=association_bugs)
#anova(lm.crp.bifidobacterium_adolescentis)
#summary(lm.crp.bifidobacterium_adolescentis)

#lm.crp.bifidobacterium_longum <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + bifidobacterium_longum + a_omega3_fs_dr*bifidobacterium_longum + (1 | participant), data=association_bugs)
#anova(lm.crp.bifidobacterium_longum)
#summary(lm.crp.bifidobacterium_longum)

#lm.crp.bifidobacterium_pseudocatenulatum <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + bifidobacterium_pseudocatenulatum + a_omega3_fs_dr*bifidobacterium_pseudocatenulatum + (1 | participant), data=association_bugs)
#anova(lm.crp.bifidobacterium_pseudocatenulatum)
#summary(lm.crp.bifidobacterium_pseudocatenulatum)

#lm.crp.bifidobacterium_bifidum <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + bifidobacterium_bifidum + a_omega3_fs_dr*bifidobacterium_bifidum + (1 | participant), data=association_bugs)
#anova(lm.crp.bifidobacterium_bifidum)
#summary(lm.crp.bifidobacterium_bifidum)

#model adjusting for other covariates, p-int=0.46 ########
lm.crp.prevotella_copri <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + prevotella_copri2c + a_omega3_fs_dr*prevotella_copri2c + calor_fs_dr + abx + (1 | participant), data=association_bugs)
anova(lm.crp.prevotella_copri)
summary(lm.crp.prevotella_copri)

#further adjusting for bristol score and physical activity, p-int=0.30 #########
lm.crp.prevotella_copri <- lmer(formula = logcrp ~ agemlvs + a_omega3_fs_dr + prevotella_copri2c + a_omega3_fs_dr*prevotella_copri2c + calor_fs_dr + abx + bristolcat + act10 + (1 | participant), data=association_bugs)
anova(lm.crp.prevotella_copri)
summary(lm.crp.prevotella_copri)

############HDL
#model adjusting for age only, p-int=0.33 ########
lm.hdl.prevotella_copri <- lmer(formula = loghdl ~ agemlvs + a_omega3_fs_dr + prevotella_copri2c + a_omega3_fs_dr*prevotella_copri2c + (1 | participant), data=association_bugs)
anova(lm.hdl.prevotella_copri)
summary(lm.hdl.prevotella_copri)

#model adjusting for other covariates, p-int=0.26 ########
lm.hdl.prevotella_copri <- lmer(formula = loghdl ~ agemlvs + a_omega3_fs_dr + prevotella_copri2c + a_omega3_fs_dr*prevotella_copri2c + calor_fs_dr + abx + (1 | participant), data=association_bugs)
anova(lm.hdl.prevotella_copri)
summary(lm.hdl.prevotella_copri)

#further adjusting for bristol score and physical activity, p-int=0.37 #########
lm.hdl.prevotella_copri <- lmer(formula = loghdl ~ agemlvs + a_omega3_fs_dr + prevotella_copri2c + a_omega3_fs_dr*prevotella_copri2c + calor_fs_dr + abx + bristolcat + act10 + (1 | participant), data=association_bugs)
anova(lm.hdl.prevotella_copri)
summary(lm.hdl.prevotella_copri)

############TG
#model adjusting for age only, p-int=0.99 ########
lm.tg.prevotella_copri <- lmer(formula = logtg ~ agemlvs + a_omega3_fs_dr + prevotella_copri2c + a_omega3_fs_dr*prevotella_copri2c + (1 | participant), data=association_bugs)
anova(lm.tg.prevotella_copri)
summary(lm.tg.prevotella_copri)

#model adjusting for other covariates, p-int=0.97 ########
lm.tg.prevotella_copri <- lmer(formula = logtg ~ agemlvs + a_omega3_fs_dr + prevotella_copri2c + a_omega3_fs_dr*prevotella_copri2c + calor_fs_dr + abx + (1 | participant), data=association_bugs)
anova(lm.tg.prevotella_copri)
summary(lm.tg.prevotella_copri)

#further adjusting for bristol score and physical activity, p-int=0.93 #########
lm.tg.prevotella_copri <- lmer(formula = logtg ~ agemlvs + a_omega3_fs_dr + prevotella_copri2c + a_omega3_fs_dr*prevotella_copri2c + calor_fs_dr + abx + bristolcat + act10 + (1 | participant), data=association_bugs)
anova(lm.tg.prevotella_copri)
summary(lm.tg.prevotella_copri)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############### Scatter plots by P copri ####################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#scatter plot with smooth line 
association_bugs$prevotella_copri2c <- as.factor(association_bugs$prevotella_copri2c)
head(association_bugs)

ggp_by_prevotella_copri2c <- ggplot(association_bugs, aes(x=a_omega3_fs_dr, y=logcrp, color=prevotella_copri2c))+ 
  geom_point(color="black", alpha = .5, shape = 21, size = 3, stroke = 1, na.rm=TRUE) +
  geom_point(aes(color=prevotella_copri2c), size=3)+
  geom_smooth(method=glm, se=FALSE, fullrange=TRUE)+
  labs(x="Dietary omega 3 intake (g/day)", y = "Log-transformed CRP")+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_text(size=20),
        legend.text = element_text(size=20), 
        axis.title=element_text(size=20,face='bold'),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.ticks =element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_blank())+
  scale_color_discrete(name ="Presence of Prevotella copri",labels=c("No (n=588)", "Yes (n=185)"))+
  annotate("text", x = 10, y = 2, label = "P-interaction=0.46", size=8)

plot(ggp_by_prevotella_copri2c)

############HDL
ggp_by_prevotella_copri2c <- ggplot(association_bugs, aes(x=a_omega3_fs_dr, y=loghdl, color=prevotella_copri2c))+ 
  geom_point(color="black", alpha = .5, shape = 21, size = 3, stroke = 1, na.rm=TRUE) +
  geom_point(aes(color=prevotella_copri2c), size=3)+
  geom_smooth(method=glm, se=FALSE, fullrange=TRUE)+
  labs(x="Dietary omega 3 intake (g/day)", y = "Log-transformed hdl")+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_text(size=20),
        legend.text = element_text(size=20), 
        axis.title=element_text(size=20,face='bold'),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.ticks =element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_blank())+
  scale_color_discrete(name ="Presence of Prevotella copri",labels=c("No (n=588)", "Yes (n=185)"))+
  annotate("text", x = 10, y = 2, label = "P-interaction=0.46", size=8)

plot(ggp_by_prevotella_copri2c)

############TG
ggp_by_prevotella_copri2c <- ggplot(association_bugs, aes(x=a_omega3_fs_dr, y=logtg, color=prevotella_copri2c))+ 
  geom_point(color="black", alpha = .5, shape = 21, size = 3, stroke = 1, na.rm=TRUE) +
  geom_point(aes(color=prevotella_copri2c), size=3)+
  geom_smooth(method=glm, se=FALSE, fullrange=TRUE)+
  labs(x="Dietary omega 3 intake (g/day)", y = "Log-transformed tg")+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_text(size=20),
        legend.text = element_text(size=20), 
        axis.title=element_text(size=20,face='bold'),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.ticks =element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_blank())+
  scale_color_discrete(name ="Presence of Prevotella copri",labels=c("No (n=588)", "Yes (n=185)"))+
  annotate("text", x = 10, y = 2, label = "P-interaction=0.46", size=8)

plot(ggp_by_prevotella_copri2c)

ggsave(filename='./ddr.crp.by.prevotella_copri2c.png', plot=ggp_by_prevotella_copri2c, width = 10, height = 5, dpi = 600) 









# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
######## Figure 1. Distribution of data and PcoA on species-level taxonomy#######
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
######### read in data########
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# read in metadata (307)
z_metadata_all <- read.table('./noGit/maaslin2/pcl/z_metadata_all.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(z_metadata_all) <- z_metadata_all$sample
z_metadata_all <- z_metadata_all[,-1]
meta_dis <- as.data.frame(t(z_metadata_all))

# read in species
species_all <- read.table('./noGit/maaslin2/pcl/species_all.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(species_all) <- species_all$sample
species_all <- species_all[,-1]
species_all <- as.data.frame(t(species_all))

# read in metadata (925)
mlvs_metadata_925 <- read.table('./noGit/maaslin2/pcl/mlvs_metadata_925.pcl', header=TRUE, sep='\t', check.names=TRUE, quote ="")
rownames(mlvs_metadata_925) <- mlvs_metadata_925 $sample
mlvs_metadata_925  <- mlvs_metadata_925 [,-1]
mlvs_metadata_925  <- as.data.frame(t(mlvs_metadata_925 ))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
####### Fig 1B. Distribution of DDR omega3 and CRP#######
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
graphics.off()

####distribution of ddr omega3, count########## Why week1?
meta_dis$omega310v
meta_dis$omega610v
meta_dis$omega3_avg

omega3_avg_distribution<-ggplot(meta_dis, aes(x=omega3_avg))+ #a_omega3_fs_dr_w1
  geom_histogram(color="black", fill="grey", binwidth=0.5)+
  labs(x="omega3_avg (g/day)", y = "Count")+
  xlim(0,5)+ #16
  theme(axis.title=element_text(size=50,face = 'bold'),
        axis.text.x = element_text(size=50),
        axis.text.y = element_text(size=50),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid')
  )
plot(omega3_avg_distribution)
ggsave(filename='./noGit/result/figures/distribution/omega3_avg_histogram_count.png', 
       plot=omega3_avg_distribution, width = 12, height = 5, dpi = 600)

omega610v_distribution<-ggplot(meta_dis, aes(x=omega610v))+ #a_omega3_fs_dr_w1
  geom_histogram(color="black", fill="grey", binwidth=2)+
  labs(x="omega610v (g/day)", y = "Count")+
  #xlim(0,5)+ #16
  theme(axis.title=element_text(size=50,face = 'bold'),
        axis.text.x = element_text(size=50),
        axis.text.y = element_text(size=50),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid')
  )
plot(omega610v_distribution)
ggsave(filename='./noGit/result/figures/distribution/omega610v_histogram_count.png', 
       plot=omega610v_distribution, width = 12, height = 5, dpi = 600)

####distribution of logcrp, count##########
summary(meta_dis$logcrp_w1)
logcrp_distribution<-ggplot(meta_dis, aes(x=logcrp_w1))+
  geom_histogram(color="black", fill="grey", binwidth=0.5)+
  labs(x="log C-reactive protein", y = "Count")+
  xlim(-2,3)+
  theme(axis.title=element_text(size=50, face = 'bold'),
        axis.text.x = element_text(size=50),
        axis.text.y = element_text(size=50),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid')
  )
plot(logcrp_distribution)
ggsave(filename='./noGit/result/figures/distribution/logcrp_histogram_count.png', plot=logcrp_distribution, width = 12, height = 5, dpi = 600)

####distribution of loghdl, count##########
meta_dis$logtg_w1

summary(meta_dis$loghdl_w1)
logcrp_distribution<-ggplot(meta_dis, aes(x=loghdl_w1))+
  geom_histogram(color="black", fill="grey", binwidth=0.25)+
  labs(x="loghdl_w1", y = "Count")+
  #xlim(0,10)+
  theme(axis.title=element_text(size=50, face = 'bold'),
        axis.text.x = element_text(size=50),
        axis.text.y = element_text(size=50),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid')
  )
plot(logcrp_distribution)
ggsave(filename='./noGit/result/figures/distribution/loghdl_histogram_count.png', plot=logcrp_distribution, width = 12, height = 5, dpi = 600)


####distribution of logtg_w1, count##########
meta_dis$logtg_w1

summary(meta_dis$logtg_w1)
logcrp_distribution<-ggplot(meta_dis, aes(x=logtg_w1))+
  geom_histogram(color="black", fill="grey", binwidth=0.25)+
  labs(x="logtg_w1", y = "Count")+
  #xlim(0,10)+
  theme(axis.title=element_text(size=50, face = 'bold'),
        axis.text.x = element_text(size=50),
        axis.text.y = element_text(size=50),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid')
  )
plot(logcrp_distribution)
ggsave(filename='./noGit/result/figures/distribution/logtg_histogram_count.png', plot=logcrp_distribution, width = 12, height = 5, dpi = 600)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
######### Fig 1C. Correlation between omega3 and age/BMI##########
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

####scatter plots of ddr omega3 and age####
ddromega3_age <- ggplot(
  data=meta_dis, aes_string(x='agemlvs', y='a_omega3_fs_dr_w1')) +
  geom_point( aes(), fill = '#C92D39', color='black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
  stat_smooth(method = "glm", color ='black')+ 
  guides(alpha='none')+labs("")+
  xlab('Age (years)') +  ylab("Dietary omega3 intake (g/day)")  +
  nature_theme +
  guides(legend.position=NULL)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(size=30, face = 'bold'),
        axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=30),
        axis.ticks =element_blank(),
        plot.title = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid')
  )+
  annotate("text", x = 73, y=15, label = "spearman r=-0.02", size=15)
plot(ddromega3_age)
ggsave(filename='./noGit/result/figures/distribution/ddromega3_age.png', plot=ddromega3_age, width = 10, height = 10, dpi = 600) 

cor.test(meta_dis$a_omega3_fs_dr_w1, meta_dis$agemlvs,
         alternative = c("two.sided"),
         method = c("spearman"), exact=F)



####scatter plots of ddr fiber and BMI####
ddrfiber_bmi <- ggplot(
  data=meta_dis, aes_string(x='bmi12', y='a_omega3_fs_dr_w1')) +
  geom_point( aes(), fill = '#C92D39', color='black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
  stat_smooth(method = "glm", color ='black')+ 
  guides(alpha='none')+labs("")+
  xlab(expression(paste("Body mass index (kg/m"^"2"*")"))) + 
  ylab("Dietary omega3 intake (g/day)")  +
  nature_theme +
  guides(legend.position=NULL)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title=element_text(size=30, face = 'bold'),
        axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=30),
        axis.ticks =element_blank(),
        plot.title = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid')
  )+
  annotate("text", x = 30, y=15, label = "spearman r=-0.24", size=15)
plot(ddrfiber_bmi)
ggsave(filename='./noGit/result/figures/distribution/ddr_omega3_bmi.png', plot=ddrfiber_bmi, width = 10, height = 10, dpi = 600) 

cor.test(meta_dis$a_omega3_fs_dr_w1, meta_dis$bmi12,
         alternative = c("two.sided"),
         method = c("spearman"), exact=F)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
######## Fig 1D. stacked barplots of fiber sources###########
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



#Exposure Data

z_mlvs_exposure_new <- read.table("noGit/mlvs_exposure_for_jorick.txt", 
                                  header=TRUE, sep='\t', check.names=TRUE, quote ="")
#z_mlvs_exposure_new <- merge(totom,z_mlvs_exposure_new,by="id")

z_idkey <- read.csv('noGit/idkey.csv', header=TRUE, sep = ",", stringsAsFactors = FALSE)
rownames(z_idkey) <- z_idkey$id
z_idkey <- z_idkey[order(z_idkey$id), ] 

z_mlvs_exposure <- reassignKey(z_mlvs_exposure_new)

#gsub search and replace in files -> ~~_w1avg colnames changed to ~~_w1, etc.
colnames(z_mlvs_exposure) <- gsub("w1avg","w1", colnames(z_mlvs_exposure)) 
colnames(z_mlvs_exposure) <- gsub("w2avg","w2", colnames(z_mlvs_exposure))
colnames(z_mlvs_exposure) <- gsub("plasma1","w1", colnames(z_mlvs_exposure))
colnames(z_mlvs_exposure) <- gsub("plasma2","w2", colnames(z_mlvs_exposure))

z_mlvs_agebmi <- z_mlvs_exposure [ , colnames(z_mlvs_exposure) %in% c('agemlvs', 'bmi12','act10','id')] 
#id is actually rownames


keep_metadata_cum <- c('ala10v','epa10v','dha10v','dpa10v','calor10v', 'trans10v', 'omega610v', 'omega310v', 'omega3_noala10v')
z_metadata_cum <- z_mlvs_exposure[, keep_metadata_cum]
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

barplot_cum <- meta_dis[,c('ala_avg','epa_avg','dha_avg','dpa_avg')]
barplot_cum <- barplot_cum %>% 
  #mutate(noala = ala_avg - epa_avg - dha_avg - dpa_avg)+
  mutate(total = ala_avg + epa_avg + dha_avg + dpa_avg)
rownames(barplot_cum) <- rownames(z_metadata_cum)
colnames(barplot_cum) <- c('ala','epa','dha','dpa', 'total')

# sort by aofib10v
barplot_cum_sort <- barplot_cum[order(-barplot_cum$total),]
order(-barplot_cum$total)

# keep only fiber subsets
barplot_cum_sort_stratified <- barplot_cum_sort[,1:4]

# transform to relative abundances
#barplot_cum_sort_stratified_a <- sweep(barplot_cum_sort_stratified, 1, rowSums(barplot_cum_sort_stratified), `/`)

# transform to samples in columns and metadata in rows
t_barplot_cum_sort_stratified <- as.data.frame(t(barplot_cum_sort_stratified))


t_barplot_cum_sort_stratified$sources <- row.names(t_barplot_cum_sort_stratified)
dim(t_barplot_cum_sort_stratified) # 4 308
# melt
t_barplot_cum_sort_stratified_melt <- melt(t_barplot_cum_sort_stratified, id.vars = 'sources')
dim(t_barplot_cum_sort_stratified_melt) # 1228 3


top30_s_colors <- scale_fill_manual(values = c("#3288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D43", "#E7298A", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142", "#D53E4F", "#FDAE61", "#FEE08B", "#333333", "#66C2A5", "#3288BD", "#1B9E77", "#D95F02", "#7570B3", "#E6AB02", "#A6761D", "#CCCCCC"))

bar_cum<-ggplot(data=t_barplot_cum_sort_stratified_melt, aes(x=variable, y=value, fill=sources)) +
  geom_bar(stat="identity")+
  theme_bw(base_size = 30) + 
  theme(plot.title = element_text(size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title=element_text(face = 'bold'),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size=30),
        axis.ticks.x=element_blank()) + 
  top30_s_colors + 
  xlab("Participant") + 
  ylab("Dietary omega3 intake (g/day)") 
#scale_y_continuous(labels = percent_format(), expand = c(0, 0)) 
print(bar_cum)

ggsave(filename='./noGit/result/figures/distribution/omega3_types.png', 
       plot=bar_cum, width = 15, height = 6, dpi = 600) 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Fig 1E. PCoA of species-level taxonomy (decorated with fruit fiber and CRP) ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
meta_pcoa <- mlvs_metadata_925[ , colnames(mlvs_metadata_925) %in% c('a_trn07_fo_dr','logcrp')] #a_omega3_fs_dr

# create quartiles for frtaf10v
quantile(meta_pcoa$a_trn07_fo_dr, na.rm=TRUE)
# 0%         25%         50%         75%        100% 
# 0.04428571  3.19000000  4.43000000  6.05571429 15.22000000 

#meta_pcoa <- meta_pcoa %>%
  #mutate(frtaf10vq = case_when(
    #0.789 < a_omega3_fs_dr & a_omega3_fs_dr <= 1.67 ~ 'Q1',
    #1.67 < a_omega3_fs_dr & a_omega3_fs_dr<= 2.21 ~ 'Q2',
    #2.21 < a_omega3_fs_dr & a_omega3_fs_dr<= 2.78 ~ 'Q3',
    #2.78 < a_omega3_fs_dr ~ 'Q4'))

meta_pcoa <- meta_pcoa %>%
  mutate(frtaf10vq = case_when(
    0.71 < a_trn07_fo_dr & a_trn07_fo_dr <= 1.30 ~ 'Q1',
    1.30 < a_trn07_fo_dr & a_trn07_fo_dr<= 1.80 ~ 'Q2',
    1.80 < a_trn07_fo_dr & a_trn07_fo_dr<= 2.41 ~ 'Q3',
    2.41 < a_trn07_fo_dr ~ 'Q4'))


table(meta_pcoa$frtaf10vq)
meta_pcoa$frtaf10vq<- factor(meta_pcoa$frtaf10vq, 
                             levels = c('Q1', 'Q2', 'Q3', 'Q4'))

# create quartiles for logcrp
quantile(meta_pcoa$logcrp, na.rm=TRUE)
# 0%        25%        50%        75%       100% 
# -1.6094379 -0.9162907 -0.1053605  0.6418539  2.9704145 

meta_pcoa <- meta_pcoa %>%
  mutate(logcrpq = case_when(
    -2.0 <= logcrp & logcrp <= -0.9162907 ~ 'Q1',
    -0.9162907 < logcrp & logcrp<= -0.1053605 ~ 'Q2',
    -0.1053605 < logcrp & logcrp<= 0.6418539 ~ 'Q3',
    0.6418539 < logcrp ~ 'Q4'))
table(meta_pcoa$logcrpq)

library(tidyverse)
meta_pcoa$logcrpq <- fct_explicit_na(meta_pcoa$logcrpq, na_level = "Missing")

rownames(meta_pcoa) <- rownames(mlvs_metadata_925)  



bugs_pcoa <- species_all
#make sure samples are in rows
bugs_pcoa <- capscale(bugs_pcoa ~ 1, distance = 'bray') #Question? old pcoa used asin sqrt transformed bug data?
variance.bugs <- head(eigenvals(bugs_pcoa)/sum(eigenvals(bugs_pcoa)))
x_variance.bugs <- as.integer(variance.bugs[1]*100)
y_variance.bugs <- as.integer(variance.bugs[2]*100)

# sample scores
ord.scores <- as.data.frame(scores(bugs_pcoa, display = "sites"))

# bug scores
ord.bug.scores <- as.data.frame(scores(bugs_pcoa, display = "species") )


# append metadata to scores
ord.scores.meta <- merge(ord.scores, meta_pcoa, by = 'row.names')

# set first column as rownames and remove it
rownames(ord.scores.meta) <- ord.scores.meta[ ,1]
ord.scores.meta[ ,1] <- NULL

# remove logCRP = NA
ord.scores.meta <- ord.scores.meta[complete.cases(ord.scores.meta), ]

# number of samples
n_samples <- nrow(ord.scores.meta)


library(viridis)
ordcolors <- scale_colour_gradientn(colours = viridis(10), limits=c(0,4))

pcoa.plot <- 
  ggplot( ord.scores.meta, aes(MDS1, MDS2) ) + 
  geom_point(size=3, alpha = 0.75, aes( color = a_trn07_fo_dr)) + #change name
  #shape = factor(logcrpq) or factor(a_omega3_fs_der)
  #scale_shape_manual(values=c(3, 16, 17, 15), name = "CRP quartiles")+
  labs(color = "Short Term Trans Fat (g/day)")+ #change name
  theme_classic() + 
  #theme_bw(base_size = 30) +
  theme(axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid'),
        axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_text( size = 30, face = 'bold'), 
        legend.text = element_text(size=20), legend.title = element_text(size=30),
        legend.position = 'right', legend.direction = 'vertical', legend.box = 'vertical',
        plot.title = element_text( size = 30, face = 'bold')
  ) +
  guides( color = guide_legend(nrow=2, byrow = TRUE) ) +
  xlab(paste0("PCo1 (",x_variance.bugs,"%)")) + ylab(paste0("PCo2 (",y_variance.bugs,"%)"))

pcoa.plot <- pcoa.plot +   ordcolors + guides(colour="colourbar")

print(pcoa.plot)

ggsave(filename='./noGit/result/figures/distribution/pcoa_bray_plot.png', plot=pcoa.plot, width = 15, height = 6, dpi = 600) 


































# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### extra code ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#met <- read.csv('/Users/jorickbater/Downloads/mlvs_exposure_wenjie.csv', header=T, row.names = 1)
#met$id = rownames(met)

#z_metadata_all$choc_avg <- (z_metadata_all$dchocffq1+z_metadata_all$mchocffq1+z_metadata_all$dchocffq2+z_metadata_all$mchocffq2)/2
#z_metadata_all$choc_avg[is.na(z_metadata_all$choc_avg)] <- 0

#Blood Fatty acids
#fa <- readxl::read_xlsx("/Users/jorickbater/Desktop/Mingyang Stuff/fattyacids.xlsx")
#fa <- as.data.frame(fa)
#fa$ID <- as.numeric(fa$ID)
#fa <- fa[,3:47]
#dat <- merge(fa,tot,by="ID")

#gut <-read_tsv("metaphlan_taxonomic_profiles.tsv")
#gut <- t(gut)
#cols <- gut[1,]
#colnames(gut) <- cols
#gut = gut[-1,]
#gut <- gut[ , grepl( "s__" , colnames(gut) ) ]
#colnames(gut) <- gsub('.*s__', '', colnames(gut))

#gut <- as.data.frame(gut)
#gut <- setDT(gut, keep.rownames = TRUE)[]
#gut$rn <- gsub("([0-9]+)_.*", "\\1", gut$rn)

#gut <- sapply(gut[,1:544], as.numeric)
#gut <- as.data.frame(gut)
#gut$rn <- as.factor(gut$rn)
#gut <- aggregate(.~rn, gut, mean)
#gut <- rename(gut,ID=rn)

#rownames(gut) = gut$ID

#z_mlvs_exposure <- read.csv('/Users/jorickbater/Desktop/Mingyang Stuff/Data/mlvs_exposure_wenjie.csv', header=T, row.names = 1)
#z_mlvs_exposure$id = rownames(z_mlvs_exposure)
#z_mlvs_exposure$id <- as.factor(z_mlvs_exposure$id)
#z_mlvs_exposure <- merge(totom,z_mlvs_exposure,by="id")

#z_metadata_all$participant <- rownames(z_metadata_all)
#met <- z_metadata_all

#Merge ID
#ID <- read_sas("/Users/jorickbater/Desktop/Mingyang Stuff/Data/mlvs_bld_phase12.sas7bdat")
#ID = ID[,c("HarvardID","id")]

#tot <- merge(met,ID, by="id")
#tot$id <- NULL
#tot <- rename(tot,id=HarvardID)

#Gut Taxa
#tax_rpk_name <-   read.table( "./bugs_dna_929_unFilt.tsv",
#                              sep = '\t',    header = TRUE,    check.names = FALSE,    na.strings = c("", "NA"))
#dim(tax_rpk_name)
#tax_rpk_name<-tax_rpk_name %>%
#  separate(Sample, c("kingdom",       "phylum",        "class" ,        "order",         "family",        "genus" ,        "species" ,      "strain"), 
#           sep = '\\|', remove = TRUE)

# only keep species-level features
#tax_rpk_species <- subset(tax_rpk_name,!is.na(species) & is.na(strain))
#rownames(tax_rpk_species)<-tax_rpk_species$species
#tax_rpk_species<-tax_rpk_species[,-c(1:8)]

#all_species <- as.data.frame(tax_rpk_species)
#all_species <- all_species[ , !colnames(all_species) %in% c(grep("15074981", colnames(all_species), value = T))]
#dim(all_species)

#dim(all_species[apply(all_species, 1, function(bug) sum(bug >= 0.01) >= 0.1*ncol(all_species)), ])
#all_species_filt <- all_species[ apply(all_species, 1, function(bug) sum(bug >= 0.01) >= 0.1*ncol(all_species)), ]
#species_list <- rownames(all_species_filt)
#species_list.format <- tolower(gsub(pattern = "s__", "", species_list))

#z_idkey <- read_sas("/Users/jorickbater/Desktop/Mingyang Stuff/Data/idlinks.sas7bdat")
#z_idkey <- z_idkey[2:909,]
#rownames(z_idkey) <- z_idkey$id
#z_idkey <- z_idkey[order(z_idkey$id), ] 
#z_idkey <- rename(z_idkey,ID1=aliasid8digits)



#med_id<-subset(tot, select=id)
#ttax_rpk <- as.data.frame(t(tax_rpk_species))
#ttax_rpk$id<-rownames(ttax_rpk)
#ttax_rpk$id <- gsub("([0-9]+)_.*", "\\1", ttax_rpk$id)

#ttax_rpk<-inner_join(med_id, ttax_rpk, by="id")
#rownames(ttax_rpk)<-ttax_rpk$id
#ttax_rpk <- subset(ttax_rpk, select = -id)
#View(ttax_rpk)

#dat <- tot
#nam = as.vector(dat$ID)
#gut <- filter(gut, ID %in% nam)

#gut$ID <- NULL
#gut[,2:544] <- gut[,2:544]/100

#rownames(dat) = dat$ID
#dat$ID = as.factor(dat$ID)

#Make Phyloseq Object
#OTU = t(gut)
#OTU = otu_table(OTU, taxa_are_rows = TRUE)
#MET = sample_data(dat, errorIfNULL=FALSE)
#phy <- phyloseq(OTU, MET)

#Filter
#library(metagMisc)
#phy <- phyloseq_filter_prevalence(phy, prev.trh = 0.10, abund.trh = 0.01, threshold_condition = "AND") 
#gut = t(phy@otu_table@.Data)
#dim(gut)

#Alpha Diversity
#plot_richness(phy, x="P3", measures = c("Shannon", "InvSimpson")) + geom_boxplot()
#plot_richness(phy, x="P10", measures = c("Shannon", "InvSimpson")) + geom_boxplot()
#plot_richness(phy, x="P13", measures = c("Shannon", "InvSimpson")) + geom_boxplot()
#plot_richness(phy, x="P14", measures = c("Shannon", "InvSimpson")) + geom_boxplot()
#plot_richness(phy, x="crp_plasma1", measures = c("Shannon", "InvSimpson")) + geom_boxplot()
#plot_richness(phy, x="omega_3_epa_ffq1", measures = c("Shannon", "InvSimpson")) + geom_boxplot()

#Beta Diversity
#source("miseqR.R")
#set.seed(1)
#dat_scale <- phy %>%
#  scale_reads(round = "round") 
#dat_bray <- phyloseq::distance(dat_scale, method = "bray")
#ordination = ordinate(phy, method="PCoA", distance=dat_bray)
#plot_ordination(dat, ordination, color="SampleType") + theme(aspect.ratio=1)

#PERMANOVA
#library(microbiome)
#library(vegan)
#rel <- microbiome::transform(phy, "compositional")

#Beta Diversity

#p <- plot_landscape(rel, method = "NMDS", distance = "bray", col = "crp_avg", size = 3)
#plot(p)

#otu <- abundances(phy)
#meta <- meta(phy)
#permanova <- adonis(t(otu) ~ P3 + P10 + P13 + P14 + crp_avg,
#data = dat, permutations=99, method = "bray")
#permanova <- adonis(t(otu) ~ crp_avg,
 #                   data = meta, permutations=999, method = "bray")

#permanova

#coef <- coefficients(permanova)["crp_avg",]
#top.coef <- coef[rev(order(abs(coef)))[1:20]]
#par(mar = c(3, 14, 2, 1))
#barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")

#Run Maaslin2
#library(Maaslin2)
#fit_data <- Maaslin2(
  #gut, dat, 'demo_output', transform = "none", heatmap_first_n = 100,
  #max_significance = 0.25,
  #fixed_effects = c('cod_liv_oil_avg','linseed_avg','flax_avg','omega_3_epa_avg', 'walnuts_avg', 'avocado_avg', 'leafvegavg', 'fishavg', 'agemlvs','calor_avg','abx_avg'),
  #random_effects = 'ID1',
  #normalization = 'NONE',
  #standardize = FALSE)


#FFQ Data
#ffq1 <- read_sas("/Users/jorickbater/Desktop/Mingyang Stuff/Data/mlvsffq1_excl_nts.sas7bdat")
#om1 <- ffq1[,c("id","pfa183n3c07_fs_ffq1", "pf205n3c07_fs_ffq1", "pf225n3c07_fs_ffq1", "pf226n3c07_fs_ffq1")]
#ffq2 <- read_sas("/Users/jorickbater/Desktop/Mingyang Stuff/Data/mlvsffq2_excl_nts.sas7bdat")
#om2 <- ffq2[,c("id","pfa183n3c07_fs_ffq2", "pf205n3c07_fs_ffq2", "pf225n3c07_fs_ffq2", "pf226n3c07_fs_ffq2")]
#dr1 <- read_sas("/Users/jorickbater/Desktop/Mingyang Stuff/Data/dr_nts_wk1_mean_e_adjusted.sas7bdat")
#dom1 <- dr1[,c("id","a_omega3_fs_dr_w1avg", "a_p22_5_fs_dr_w1avg")]
#dr2 <- read_sas("/Users/jorickbater/Desktop/Mingyang Stuff/Data/dr_nts_wk2_mean_e_adjusted.sas7bdat")
#dom2 <- dr2[,c("id","a_omega3_fs_dr_w2avg", "a_p22_5_fs_dr_w2avg")]

#dat_list = list(om1,om2,dom1,dom2)

#my_merge <- function(df1, df2){                               
#  merge(df1, df2, by = "id")
#}

#totom <- Reduce(my_merge,dat_list)
#totom$pfa183n3c07_fs_avg <- rowMeans(totom[,c('pfa183n3c07_fs_ffq1', 'pfa183n3c07_fs_ffq2')], na.rm=TRUE)
#totom$pf205n3c07_fs_avg <- rowMeans(totom[,c('pf205n3c07_fs_ffq1', 'pf205n3c07_fs_ffq2')], na.rm=TRUE)
#totom$pf225n3c07_fs_avg <- rowMeans(totom[,c('pf225n3c07_fs_ffq1', 'pf225n3c07_fs_ffq2')], na.rm=TRUE)
#totom$pf226n3c07_fs_avg <- rowMeans(totom[,c('pf226n3c07_fs_ffq1', 'pf226n3c07_fs_ffq2')], na.rm=TRUE)

