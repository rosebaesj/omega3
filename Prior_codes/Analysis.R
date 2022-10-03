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
mlvs_metadata_925 <- read.table('./noGit/maaslin2/pcl/mlvs_metadata_925.pcl', header=TRUE, sep='\t', check.names=FALSE, quote ="")
rownames(mlvs_metadata_925) <- mlvs_metadata_925 $sample
mlvs_metadata_925  <- mlvs_metadata_925 [,-1]
mlvs_metadata_925  <- as.data.frame(t(mlvs_metadata_925 ))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
####### +  Fig 1B. Distribution of DDR omega3 and CRP#######
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
graphics.off()

#### ++  of ddr omega3, count#########
meta_dis$omega310v
meta_dis$omega610v
meta_dis$omega3_avg
meta_dis$a_omega3_fs_dr_w1

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
ggsave(filename='./noGit/result/figures/distribution/histogram/omega3_avg_histogram_count.png', 
       plot=omega3_avg_distribution, width = 12, height = 5, dpi = 600)


omega3_avg_distribution<-ggplot(meta_dis, aes(x=omega3_noala_avg))+ #a_omega3_fs_dr_w1
  geom_histogram(color="black", fill="grey", binwidth=0.25)+
  labs(x="omega3_noala_avg (g/day)", y = "Count")+
  xlim(0,2.5)+ #16
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
ggsave(filename='./noGit/result/figures/distribution/histogram/omega3_noala_avg_histogram_count.png', 
       plot=omega3_avg_distribution, width = 12, height = 5, dpi = 600)


#a_omega3_fs_dr_w1

omega3_avg_distribution<-ggplot(meta_dis, aes(x=a_omega3_fs_dr_w1))+ #a_omega3_fs_dr_w1
  geom_histogram(color="black", fill="grey", binwidth=0.5)+
  labs(x="a_omega3_fs_dr_w1 (g/day)", y = "Count")+
  xlim(0,6)+ #16
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
ggsave(filename='./noGit/result/figures/distribution/histogram/a_omega3_fs_dr_w1_histogram_count2.png', 
       plot=omega3_avg_distribution, width = 12, height = 5, dpi = 600)


omega3_avg_distribution<-ggplot(meta_dis, aes(x=a_omega3_noala_fs_dr_w1))+ #a_omega3_fs_dr_w1
  geom_histogram(color="black", fill="grey", binwidth=0.5)+
  labs(x="a_omega3_noala_fs_dr_w1 (g/day)", y = "Count")+
  xlim(0,6)+ #16
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
ggsave(filename='./noGit/result/figures/distribution/histogram/a_omega3_noala_fs_dr_w1_histogram_count2.png', 
       plot=omega3_avg_distribution, width = 12, height = 5, dpi = 600)


omega3_avg_distribution<-ggplot(meta_dis, aes(x=a_omega3_fs_dr_w1))+ #a_omega3_fs_dr_w1
  geom_histogram(color="black", fill="grey", binwidth=0.5)+
  labs(x="a_omega3_fs_dr_w1 (g/day)", y = "Count")+
  xlim(0,6)+ #16
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
ggsave(filename='./noGit/result/figures/distribution/histogram/a_omega3_fs_dr_w1_histogram_count2.png', 
       plot=omega3_avg_distribution, width = 12, height = 5, dpi = 600)



#### ++  of ddr omega3, count#########
omega610v_distribution<-ggplot(meta_dis, aes(x=omega610v))+ #a_omega3_fs_dr_w1
  geom_histogram(color="black", fill="grey", binwidth=1.5)+
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
ggsave(filename='./noGit/result/figures/distribution/histogram/omega610v_histogram_count.png', 
       plot=omega610v_distribution, width = 12, height = 5, dpi = 600)












#### ++  of logcrp, count##########
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
ggsave(filename='./noGit/result/figures/distribution/histogram/logcrp_histogram_count.png', plot=logcrp_distribution, width = 12, height = 5, dpi = 600)


#### ++  of loghdl, count##########

summary(meta_dis$loghdl_w1)
logcrp_distribution<-ggplot(meta_dis, aes(x=loghdl_w1))+
  geom_histogram(color="black", fill="grey", binwidth=0.2)+
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
ggsave(filename='./noGit/result/figures/distribution/histogram/loghdl_histogram_count.png', plot=logcrp_distribution, width = 12, height = 5, dpi = 600)


#### ++  of logtg_w1, count##########
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
ggsave(filename='./noGit/result/figures/distribution/histogram/logtg_histogram_count.png', plot=logcrp_distribution, width = 12, height = 5, dpi = 600)




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
######### Fig 1C. Correlation between omega3 and age/BMI##########
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#### ++ ddr omega3 and age, scatter plot ####
ddromega3_age <- ggplot(
  data=meta_dis, aes_string(x='agemlvs', y='a_omega3_fs_dr_w1')) +
  geom_point( aes(), fill = '#C92D39', color='black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
  stat_smooth(method = "glm", color ='black')+ 
  guides(alpha='none')+labs("")+
  xlab('Age (years)') +  ylab("Dietary omega3 intake (g/day)")  +
  nature_theme +
  guides(legend.position=NULL)+
  ylim(0,7.5)+
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
  annotate("text", x = 73, y=7.5, label = "spearman r=-0.02", size=15)
plot(ddromega3_age)
ggsave(filename='./noGit/result/figures/distribution/scatterplot/ddromega3_age.png', plot=ddromega3_age, width = 10, height = 10, dpi = 600) 

cor.test(meta_dis$a_omega3_fs_dr_w1, meta_dis$agemlvs,
         alternative = c("two.sided"),
         method = c("spearman"), exact=F)




ddromega3_age <- ggplot(
  data=meta_dis, aes_string(x='agemlvs', y='a_omega3_noala_fs_dr_w1')) +
  geom_point( aes(), fill = '#C92D39', color='black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
  stat_smooth(method = "glm", color ='black')+ 
  guides(alpha='none')+labs("")+
  xlab('Age (years)') +  ylab("Dietary omega3 noala intake (g/day)")  +
  nature_theme +
  guides(legend.position=NULL)+
  ylim(0,2.5)+
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
  annotate("text", x = 73, y=7.5, label = "spearman r=-0.107", size=15)
plot(ddromega3_age)
ggsave(filename='./noGit/result/figures/distribution/scatterplot/ddromega3_age_noala.png', plot=ddromega3_age, width = 10, height = 10, dpi = 600) 

cor.test(meta_dis$a_omega3_noala_fs_dr_w1, meta_dis$agemlvs,
         alternative = c("two.sided"),
         method = c("spearman"), exact=F)






#### ++ ddr omega3 and BMI scatter plot ####

cor.test(meta_dis$a_omega3_fs_dr_w1, meta_dis$bmi12,
         alternative = c("two.sided"),
         method = c("spearman"), exact=F)


ddrfiber_bmi <- ggplot(
  data=meta_dis, aes_string(x='bmi12', y='a_omega3_fs_dr_w1')) +
  geom_point( aes(), fill = '#C92D39', color='black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
  stat_smooth(method = "glm", color ='black')+ 
  guides(alpha='none')+labs("")+
  xlab(expression(paste("Body mass index (kg/m"^"2"*")"))) + 
  ylab("Dietary omega3 intake (g/day)")  +
  nature_theme +
  guides(legend.position=NULL)+
  ylim(0,7.5)+
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
  annotate("text", x = 30, y=7.5, label = "spearman r=-0.11", size=15)
plot(ddrfiber_bmi)
ggsave(filename='./noGit/result/figures/distribution/scatterplot/ddr_omega3_bmi.png', plot=ddrfiber_bmi, width = 10, height = 10, dpi = 600) 








cor.test(meta_dis$a_omega3_noala_fs_dr_w1, meta_dis$bmi12,
         alternative = c("two.sided"),
         method = c("spearman"), exact=F)


ddrfiber_bmi <- ggplot(
  data=meta_dis, aes_string(x='bmi12', y='a_omega3_noala_fs_dr_w1')) +
  geom_point( aes(), fill = '#C92D39', color='black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
  stat_smooth(method = "glm", color ='black')+ 
  guides(alpha='none')+labs("")+
  xlab(expression(paste("Body mass index (kg/m"^"2"*")"))) + 
  ylab("Dietary omega3 noala intake (g/day)")  +
  nature_theme +
  guides(legend.position=NULL)+
  ylim(0,2.5)+
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
  annotate("text", x = 30, y=7.5, label = "spearman r=-0.105", size=15)
plot(ddrfiber_bmi)
ggsave(filename='./noGit/result/figures/distribution/scatterplot/ddr_omega3_bmi_noala.png', plot=ddrfiber_bmi, width = 10, height = 10, dpi = 600) 


#### ++ 3 ddr omega3 and total energy intake???? scatter plot ####
# 
# cor.test(meta_dis$a_omega3_fs_dr_w1, meta_dis$bmi12,
#          alternative = c("two.sided"),
#          method = c("spearman"), exact=F)
# 
# 
# ddrfiber_bmi <- ggplot(
#   data=meta_dis, aes_string(x='bmi12', y='a_omega3_fs_dr_w1')) +
#   geom_point( aes(), fill = '#C92D39', color='black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
#   stat_smooth(method = "glm", color ='black')+ 
#   guides(alpha='none')+labs("")+
#   xlab(expression(paste("Body mass index (kg/m"^"2"*")"))) + 
#   ylab("Dietary omega3 intake (g/day)")  +
#   nature_theme +
#   guides(legend.position=NULL)+
#   ylim(0,7.5)+
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.title=element_text(size=30, face = 'bold'),
#         axis.text.x = element_text(size=30),
#         axis.text.y = element_text(size=30),
#         axis.ticks =element_blank(),
#         plot.title = element_blank(),
#         axis.line.x = element_line(colour = 'black', size=0.75, linetype='solid'),
#         axis.line.y = element_line(colour = 'black', size=0.75, linetype='solid')
#   )+
#   annotate("text", x = 30, y=7.5, label = "spearman r=-0.11", size=15)
# plot(ddrfiber_bmi)
# ggsave(filename='./noGit/result/figures/distribution/scatterplot/ddr_omega3_bmi.png', plot=ddrfiber_bmi, width = 10, height = 10, dpi = 600) 
# 
# 



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
######## Fig 1D. stacked barplots of omega3 types###########
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

##### + intake short-term ########

## a_omega3_noala_fs_dr_w1avg
murin_omega3 <- c('a_omega3_fs_dr_w1', 'a_omega3_noala_fs_dr_w1', a_omega3_noala_fs_dr_w1)
z_murin_omega3 <- z_mlvs_exposure



keep_metadata_cum <- c('ala10v','epa10v','dha10v','dpa10v','calor10v', 'trans10v', 'omega610v', 'omega310v', 'omega3_noala10v')
z_metadata_cum <- z_mlvs_exposure[, keep_metadata_cum]
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

barplot_cum <- meta_dis[,c('ala_avg','epa_avg','dha_avg','dpa_avg')]
barplot_cum <- barplot_cum %>% 
  #mutate(noala = ala_avg - epa_avg - dha_avg - dpa_avg)+
  # mutate(total = ala_avg + epa_avg + dha_avg + dpa_avg)
  mutate(total_marine =  epa_avg + dha_avg + dpa_avg)
rownames(barplot_cum) <- rownames(z_metadata_cum)
colnames(barplot_cum) <- c('ala','epa','dha','dpa', 'total_marine')

# sort by aofib10v
barplot_cum_sort <- barplot_cum[order(-barplot_cum$total),]
order(-barplot_cum$total)

# keep only fiber subsets
barplot_cum_sort_stratified <- barplot_cum_sort[,2:4]

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


##### + plasma w1 ########
meta_dis_p <- as.data.frame(z_mlvs_exposure_plasma)

barplot_pfa <- meta_dis_p[,c('ala_pfa','epa_pfa','dha_pfa','dpa_pfa')]
barplot_pfa <- barplot_pfa %>% 
  #mutate(noala = ala_avg - epa_avg - dha_avg - dpa_avg)+
  # mutate(total = ala_avg + epa_avg + dha_avg + dpa_avg)
  mutate(total_marine =  epa_pfa + dha_pfa + dpa_pfa)
rownames(barplot_pfa) <- rownames(z_mlvs_exposure_plasma)
colnames(barplot_pfa) <- c('ala','epa','dha','dpa', 'total_marine')

# sort by aofib10v
barplot_pfa_sort <- barplot_pfa[order(-barplot_pfa$total),]
order(-barplot_pfa$total)

# keep only fiber subsets
barplot_pfa_sort_stratified <- barplot_pfa_sort[,2:4]

# transform to relative abundances
#barplot_pfa_sort_stratified_a <- sweep(barplot_pfa_sort_stratified, 1, rowSums(barplot_pfa_sort_stratified), `/`)

# transform to samples in columns and metadata in rows


t_barplot_pfa_sort_stratified <- as.data.frame(t(barplot_pfa_sort_stratified))


t_barplot_pfa_sort_stratified$sources <- row.names(t_barplot_pfa_sort_stratified)
dim(t_barplot_pfa_sort_stratified) # 4 308
# melt
t_barplot_pfa_sort_stratified_melt <- melt(t_barplot_pfa_sort_stratified, id.vars = 'sources')
dim(t_barplot_pfa_sort_stratified_melt) # 1228 3


top30_s_colors <- scale_fill_manual(values = c("#3288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D43", "#E7298A", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142", "#D53E4F", "#FDAE61", "#FEE08B", "#333333", "#66C2A5", "#3288BD", "#1B9E77", "#D95F02", "#7570B3", "#E6AB02", "#A6761D", "#CCCCCC"))

bar_pfa<-ggplot(data=t_barplot_pfa_sort_stratified_melt, aes(x=variable, y=value, fill=sources)) +
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
  ylab("Plasma omega3") 
#scale_y_continuous(labels = percent_format(), expand = c(0, 0)) 
print(bar_pfa)

ggsave(filename='./noGit/result/figures/distribution/omega3_plasma.png', 
       plot=bar_pfa, width = 15, height = 6, dpi = 600) 



##### + plasma w1 ########
meta_dis_p <- as.data.frame(z_mlvs_exposure_plasma)

barplot_pfa <- meta_dis_p[,c('ala_pfa','epa_pfa','dha_pfa','dpa_pfa')]
barplot_pfa <- barplot_pfa %>% 
  #mutate(noala = ala_avg - epa_avg - dha_avg - dpa_avg)+
  # mutate(total = ala_avg + epa_avg + dha_avg + dpa_avg)
  mutate(total_marine =  ala_pfa + epa_pfa + dha_pfa + dpa_pfa)
rownames(barplot_pfa) <- rownames(z_mlvs_exposure_plasma)
colnames(barplot_pfa) <- c('ala','epa','dha','dpa', 'total_marine')

# sort by aofib10v
barplot_pfa_sort <- barplot_pfa[order(-barplot_pfa$total),]
order(-barplot_pfa$total)

# keep only fiber subsets
barplot_pfa_sort_stratified <- barplot_pfa_sort[,1:4]

# transform to relative abundances
#barplot_pfa_sort_stratified_a <- sweep(barplot_pfa_sort_stratified, 1, rowSums(barplot_pfa_sort_stratified), `/`)

# transform to samples in columns and metadata in rows


t_barplot_pfa_sort_stratified <- as.data.frame(t(barplot_pfa_sort_stratified))


t_barplot_pfa_sort_stratified$sources <- row.names(t_barplot_pfa_sort_stratified)
dim(t_barplot_pfa_sort_stratified) # 4 308
# melt
t_barplot_pfa_sort_stratified_melt <- melt(t_barplot_pfa_sort_stratified, id.vars = 'sources')
dim(t_barplot_pfa_sort_stratified_melt) # 1228 3


top30_s_colors <- scale_fill_manual(values = c("#3288BD", "#5E4FA2", "#66A61E", "#FFED6F", "#FF0000", "#F46D43", "#E7298A", "#00008B", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5", "#9E0142", "#D53E4F", "#FDAE61", "#FEE08B", "#333333", "#66C2A5", "#3288BD", "#1B9E77", "#D95F02", "#7570B3", "#E6AB02", "#A6761D", "#CCCCCC"))

bar_pfa<-ggplot(data=t_barplot_pfa_sort_stratified_melt, aes(x=variable, y=value, fill=sources)) +
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
  ylab("Plasma omega3") 
#scale_y_continuous(labels = percent_format(), expand = c(0, 0)) 
print(bar_pfa)

ggsave(filename='./noGit/result/figures/distribution/omega3_plasma_withala.png', 
       plot=bar_pfa, width = 15, height = 6, dpi = 600) 






# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#### Fig 1E. PCoA of species-level taxonomy ####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
library(viridis)






##### + total trans-fatty acids #####
#### dot by fecal sample which is ~4 times participant
meta_pcoa <- mlvs_metadata_925[ , colnames(mlvs_metadata_925) %in% c('a_trn07_fo_dr','logcrp')] #a_omega3_fs_dr
# create quartiles for frtaf10v
quantile(meta_pcoa$a_trn07_fo_dr, na.rm=TRUE)
# 0%         25%         50%         75%        100% 
# 0.04428571  3.19000000  4.43000000  6.05571429 15.22000000 


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
  scale_shape_manual(values=c(3, 16, 17, 15), name = "CRP quartiles")+
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

ggsave(filename='./noGit/result/figures/distribution/pcoa/short_transfat.png', plot=pcoa.plot, width = 15, height = 6, dpi = 600) 



##### + total marine omega3 #####
#### dot by fecal sample which is ~4 times participant
meta_pcoa <- mlvs_metadata_925[ , colnames(mlvs_metadata_925) %in% c('a_omega3_noala_fs_dr','a_omega3_fs_dr','logcrp')] #a_omega3_fs_dr
# create quartiles

quantile(meta_pcoa$a_omega3_noala_fs_dr, na.rm=TRUE)
#          0%         25%         50%         75%        100% 
#0.002604295 0.144673841 0.412159454 0.810939233 6.535326310 

meta_pcoa <- meta_pcoa %>%
  mutate(a_omega3_noala_fs_dr_Q = case_when(
    a_omega3_noala_fs_dr < 0.144673841 ~ 'Q1',
    0.144673841 <= a_omega3_noala_fs_dr & a_omega3_noala_fs_dr< 0.412159454 ~ 'Q2',
    0.412159454 <= a_omega3_noala_fs_dr & a_omega3_noala_fs_dr< 0.810939233 ~ 'Q3',
    0.810939233 <= a_omega3_noala_fs_dr ~ 'Q4'))


table(meta_pcoa$a_omega3_noala_fs_dr_Q)
meta_pcoa$a_omega3_noala_fs_dr_Q<- factor(meta_pcoa$a_omega3_noala_fs_dr_Q, 
                             levels = c('Q1', 'Q2', 'Q3', 'Q4'))


quantile(meta_pcoa$a_omega3_fs_dr, na.rm=TRUE)
#       0%        25%        50%        75%       100% 
# 0.7894288  1.6667867  2.2143668  2.7808156 12.7327176 

meta_pcoa <- meta_pcoa %>%
  mutate(a_omega3_fs_dr_Q = case_when(
    a_omega3_fs_dr < 1.6667867 ~ 'Q1',
    1.6667867 <= a_omega3_fs_dr & a_omega3_fs_dr< 2.2143668 ~ 'Q2',
    2.2143668 <= a_omega3_fs_dr & a_omega3_fs_dr< 2.7808156 ~ 'Q3',
    2.7808156 <= a_omega3_fs_dr ~ 'Q4'))
table(meta_pcoa$a_omega3_fs_dr_Q)
meta_pcoa$a_omega3_fs_dr_Q<- factor(meta_pcoa$a_omega3_fs_dr_Q, 
                                          levels = c('Q1', 'Q2', 'Q3', 'Q4'))

# create quartiles for logcrp
quantile(meta_pcoa$logcrp, na.rm=TRUE)
# 0%        25%        50%        75%       100% 
# -1.6094379 -0.9162907 -0.1053605  0.6418539  2.9704145 

meta_pcoa <- meta_pcoa %>%
  mutate(logcrp_Q = case_when(
    logcrp <= -0.9162907 ~ 'Q1',
    -0.9162907 < logcrp & logcrp<= -0.1053605 ~ 'Q2',
    -0.1053605 < logcrp & logcrp<= 0.6418539 ~ 'Q3',
    0.6418539 < logcrp ~ 'Q4'))
table(meta_pcoa$logcrp_Q)

meta_pcoa$logcrp_Q<- factor(meta_pcoa$logcrp_Q, 
                                    levels = c('Q1', 'Q2', 'Q3', 'Q4'))





library(tidyverse)
#meta_pcoa$logcrp_Q <- fct_explicit_na(meta_pcoa$logcrp_Q, na_level = "Missing")

rownames(meta_pcoa) <- rownames(mlvs_metadata_925)  



bugs_pcoa <- species_all
#make sure samples are in rows
bugs_pcoa <- capscale(bugs_pcoa ~ 1, distance = 'bray') 
#Question? old pcoa used asin sqrt transformed bug data?
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




ggplot( ord.scores.meta, aes(MDS1, MDS2) ) + 
  geom_point(size=3, alpha = 0.75, aes( color = a_omega3_noala_fs_dr, shape = logcrp_Q))+#, shape = logcrp_Q)) + #change name
  #shape = factor(logcrp_Q) +#or factor(a_omega3_fs_der)
  scale_shape_manual(values=c(16, 17, 15, 18), name = "CRP quartiles")+
  labs(color = "Short Term Omega3 no ala (g/day)")+ #change name
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
  scale_colour_gradientn(colours = viridis(10), limits=c(0,1))+
  guides(colour="colourbar")+
  xlab(paste0("PCo1 (",x_variance.bugs,"%)")) + ylab(paste0("PCo2 (",y_variance.bugs,"%)"))


ggsave(filename='./noGit/result/figures/distribution/pcoa/short_woala.png', 
       width = 13, height = 6, dpi = 600) 




ggplot( ord.scores.meta, aes(MDS1, MDS2) ) + 
  geom_point(size=3, alpha = 0.75, aes( color = a_omega3_fs_dr, shape = logcrp_Q))+#, shape = logcrp_Q)) + #change name
  #shape = factor(logcrp_Q) +#or factor(a_omega3_fs_der)
  scale_shape_manual(values=c(16, 17, 15, 18), name = "CRP quartiles")+
  labs(color = "Short Term Omega3 with ala (g/day)")+ #change name
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
  scale_colour_gradientn(colours = viridis(10), limits=c(0,6))+
  guides(colour="colourbar")+
  xlab(paste0("PCo1 (",x_variance.bugs,"%)")) + ylab(paste0("PCo2 (",y_variance.bugs,"%)"))


ggsave(filename='./noGit/result/figures/distribution/pcoa/short_wala.png', 
       width = 13, height = 6, dpi = 600) 



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
####################### Table 1 Characteristics ################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


##### + total omega3 #####

z_mlvs_quartile <- data.frame(t(z_metadata_all))



# create quartiles
quantile(z_mlvs_quartile$a_omega3_fs_dr_w1, na.rm=TRUE)
#        0%        25%        50%        75%       100% 
# 0.5716714  1.6168791  2.1443837  2.7440490 12.7327176 

z_mlvs_quartile <- z_mlvs_quartile %>%
  mutate(a_omega3_fs_dr_w1_Q = case_when(
    a_omega3_fs_dr_w1 < 1.6168791 ~ 'Q1',
    1.6168791 <= a_omega3_fs_dr_w1 & a_omega3_fs_dr_w1< 2.1443837 ~ 'Q2',
    2.1443837 <= a_omega3_fs_dr_w1 & a_omega3_fs_dr_w1< 2.7440490 ~ 'Q3',
    2.7440490 <= a_omega3_fs_dr_w1 ~ 'Q4'))

quantile(z_mlvs_quartile$a_omega3_noala_fs_dr_w1, na.rm=TRUE)
#          0%         25%         50%         75%        100% 
#0.002604295 0.144673841 0.412159454 0.810939233 6.535326310 

z_mlvs_quartile <- z_mlvs_quartile %>%
  mutate(a_omega3_noala_fs_dr_w1_Q = case_when(
    a_omega3_noala_fs_dr_w1 < 0.144673841 ~ 'Q1',
    0.144673841 <= a_omega3_noala_fs_dr_w1 & a_omega3_noala_fs_dr_w1< 0.412159454 ~ 'Q2',
    0.412159454 <= a_omega3_noala_fs_dr_w1 & a_omega3_noala_fs_dr_w1< 0.810939233 ~ 'Q3',
    0.810939233 <= a_omega3_noala_fs_dr_w1 ~ 'Q4'))



#925 data
#         0%         25%         50%         75%        100% 
# 0.002604295 0.148045863 0.417087463 0.839926247 5.656422577 

# z_mlvs_quartile <- z_mlvs_quartile %>%
#   mutate(a_omega3_fs_dr_w1avg_Q = case_when(
#     a_omega3_fs_dr_w1avg < 0.148045863 ~ 'Q1',
#     0.148045863 <= a_omega3_fs_dr_w1avg & a_omega3_fs_dr_w1avg< 0.417087463 ~ 'Q2',
#     0.417087463 <= a_omega3_fs_dr_w1avg & a_omega3_fs_dr_w1avg< 0.839926247 ~ 'Q3',
#     0.839926247 <= a_omega3_fs_dr_w1avg ~ 'Q4'))




table(z_mlvs_quartile$a_omega3_fs_dr_w1_Q)
table(z_mlvs_quartile$a_omega3_noala_fs_dr_w1_Q)

z_mlvs_quartile$a_omega3_fs_dr_w1_Q<- factor(z_mlvs_quartile$a_omega3_fs_dr_w1_Q, 
                             levels = c('Q1', 'Q2', 'Q3', 'Q4'))
z_mlvs_quartile$a_omega3_noala_fs_dr_w1_Q<- factor(z_mlvs_quartile$a_omega3_noala_fs_dr_w1_Q, 
                                             levels = c('Q1', 'Q2', 'Q3', 'Q4'))


characteristics <- 
  z_mlvs_quartile %>% 
  group_by(a_omega3_fs_dr_w1_Q, n) %>%
  dplyr::summarise_at(funs(list( mean=mean,
                              sd=sd), na.rm=TRUE))
    
characteristics_noala <- 
  z_mlvs_quartile %>% 
  group_by(a_omega3_noala_fs_dr_w1_Q) %>%
  dplyr::summarize_at( vars(-a_omega3_fs_dr_w1_Q), list(mean=mean,
                              sd=sd), na.rm=TRUE)
    
    
t(characteristics)
write.table(t(characteristics), file = "noGit/result/characteristics_omega3_Q.tsv", sep = "\t")
write.table(t(characteristics_noala), file = "noGit/result/characteristics_noala_Q.tsv", sep = "\t")








































