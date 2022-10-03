##########
# interaction
getwd()
source("noGit/Rstart.R")
library("lme4")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
######### Figure 2. interaction ##***between DDR fiber and P. copri on logCRP****#####
# -> between 
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
############### Scatter plots by  ####################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#scatter plot with smooth line 
association_bugs$prevotella_copri2c <- as.factor(association_bugs$prevotella_copri2c)
head(association_bugs)

ggp_by_prevotella_copri2c <- ggplot(association_bugs, aes(x=a_omega3_fs_dr, y=logcrp, color=prevotella_copri2c))+ 
  geom_point(color="black", alpha = .5, shape = 21, size = 3, stroke = 1, na.rm=TRUE) +
  geom_point(aes(color="pink"), size=3)+
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
  #scale_color_discrete(name ="Presence of Prevotella copri",labels=c("No (n=588)", "Yes (n=185)"))+
  annotate("text", x = 10, y = 2, label = "P-interaction=0.46", size=8)

plot(ggp_by_prevotella_copri2c)



ggp_by_omega3_noala <- ggplot(association_bugs, aes(x=a_omega3_noala_fs_dr, y=logcrp))+ 
  geom_point(color="black", alpha = .5, shape = 21, size = 3, stroke = 1, na.rm=TRUE) +
  #geom_point(aes(), color="pink", size=3)+
  geom_smooth(method=glm, se=FALSE, fullrange=TRUE)+
  labs(x="Dietary omega 3 no ala intake (g/day)", y = "Log-transformed CRP")+
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
  #scale_color_discrete(name ="Presence of Prevotella copri",labels=c("No (n=588)", "Yes (n=185)"))+
  annotate("text", x = 5, y = 2, label = "P-interaction=0.46", size=8)

plot(ggp_by_omega3_noala)

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




