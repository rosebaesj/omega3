#########
#Figure 3






# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
####################### Fig 3A. MaAsLin2 ################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################# read in data ###################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# read in metadata (925)
mlvs_metadata_925 <- read.table('./noGit/maaslin2/pcl/mlvs_metadata_925.pcl', header=TRUE, sep='\t', 
                                check.names=FALSE, quote ="")
rownames(mlvs_metadata_925) <- mlvs_metadata_925$sample
mlvs_metadata_925  <- mlvs_metadata_925 [,-1]
mlvs_metadata_925  <- as.data.frame(t(mlvs_metadata_925 ))
rownames(mlvs_metadata_925)<-sapply(str_remove_all(rownames(mlvs_metadata_925),"X"),"[")

# read in species
species_all <- read.table('./noGit/maaslin2/pcl/species_all.pcl', header=TRUE, sep='\t', 
                          check.names=FALSE, quote ="")
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



aa <- read.table(paste('./noGit/maaslin2/pcl/a_omega3_fs_ddr.pcl', sep=''),
           header = TRUE, row.names = 1
           #  ,check.names = FALSE ##this gets rid of Xs in colnames
)




bug_file<-c('species_all')
metadata_file <- c( 'ala_ffqcum','epa_ffqcum','dha_ffqcum', 'dpa_ffqcum', 'trans_ffqcum', 'omega6_ffqcum', 'omega3_ffqcum', 'omega3_noala_ffqcum',
                    'ala_avg_fs_ffq', 'epa_avg_fs_ffq', 'dha_avg_fs_ffq', 'dpa_avg_fs_ffq', 'trans_avg_fs_ffq', 'omega6_avg_fs_ffq', 'omega3_avg_fs_ffq', 'omega3_noala_avg_fs_ffq', 'a_omega3_fs_ddr', 'a_omega3_noala_fs_ddr', 'a_trans_fs_ddr', 'a_ala_fs_ddr', 'a_epa_fs_ddr', 'a_dha_fs_ddr', 'a_dpa_fs_ddr',
                    'biomarker_logcrp', 'biomarker_logcrp_complete', #'a_omega3_fs_ddr_logcrp', 
                    'biomarker_loghdl', 'biomarker_loghdl_complete', #'a_omega3_fs_ddr_loghdl', 
                    'biomarker_logtg', 'biomarker_logtg_complete'#, #'a_omega3_fs_ddr_logtg', 
                    # 'a_trans_fs_ddr_logcrp', 'a_trans_fs_ddr_loghdl', 'a_trans_fs_ddr_logtg',
                    # 'a_ala_fs_ddr_logcrp', 'a_ala_fs_ddr_loghdl', 'a_ala_fs_ddr_logtg',
                    # 'a_epa_fs_ddr_logcrp', 'a_epa_fs_ddr_loghdl', 'a_epa_fs_ddr_logtg',
                    # 'a_dha_fs_ddr_logcrp', 'a_dha_fs_ddr_loghdl', 'a_dha_fs_ddr_logtg',
                    # 'a_dpa_fs_ddr_logcrp', 'a_dpa_fs_ddr_loghdl', 'a_dpa_fs_ddr_logtg'
                    )

pfa_file <- c('ala_pfa', 'epa_pfa', 'dpa_pfa', 'dha_pfa')

input_metadata <- read.table(paste('./noGit/maaslin2/pcl/biomarker_logcrp.pcl', sep=''),
                             header = TRUE, row.names = 1
                             #  ,check.names = FALSE ##this gets rid of Xs in colnames
)

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

for (a in 1:length(pfa_file)){
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
             output           = paste('./noGit/maaslin2/plasmafa', metadata_file[a], '_', bug_file[b], '/', sep=''),
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

colnames(ddr_plot[,c(1:41)])
features <- ddr_plot[,c(1:41)]
microbiomes <- ddr_plot[,c(42:ncol(ddr_plot))]

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


adjust <- ddr_plot[,c('agemlvs', 'abx', 'calor_fs_dr','participant')]
adjuststool <- ddr_plot[,c('agemlvs', 'abx', 'calor_fs_dr','participant')]
a <- cbind(adjust, a_f205_fs_dr = ddr_plot[,c('a_f205_fs_dr')])

murine_omega3 <- ddr_plot[,c('participant',
  'a_f205_fs_dr', 'a_f226_fs_dr', 'a_p22_5_fs_dr','a_omega3_noala_fs_dr')]


Maaslin2(input_data = microbiomes, 
         input_metadata = murine_omega3,
         output           = './noGit/maaslin2/a_murine_omega3_log/',
         normalization    = 'NONE', 
         standardize      = 'FALSE',
         transform        = "LOG",#'AST', 
         analysis_method  = 'LM', 
         random_effects   = 'participant',
         min_abundance    = 0, 
         min_prevalence   = 0,
         cores            = 1)

Maaslin2(input_data = microbiomes, 
         input_metadata = murine_omega3,
         output           = './noGit/maaslin2/a_murine_omega3_CPLM/',
         normalization    = 'NONE', 
         standardize      = 'FALSE',
         transform        = 'NONE', #CPLM can;t AST because AST makes negative values
         analysis_method  = 'CPLM', #instead of LM because many values are close to zero
         random_effects   = 'participant',
         min_abundance    = 0, 
         min_prevalence   = 0,
         cores            = 1)

Maaslin2(input_data = microbiomes, 
         input_metadata = cbind(adjust, a_f205_fs_dr = ddr_plot[,c('a_f205_fs_dr')]),
         output           = './noGit/maaslin2/EPA_a_f205_fs_dr/',
         normalization    = 'NONE', 
         standardize      = 'FALSE',
         transform        = 'AST', 
         analysis_method  = 'LM', 
         random_effects   = 'participant',
         min_abundance    = 0, 
         min_prevalence   = 0,
         cores            = 1)

Maaslin2(input_data = microbiomes, 
         input_metadata = ddr_plot[,c('participant','a_f226_fs_dr')],
         output           = './noGit/maaslin2/DHA_a_f226_fs_dr/',
         normalization    = 'NONE', 
         standardize      = 'FALSE',
         transform        = 'AST', 
         analysis_method  = 'LM', 
         random_effects   = 'participant',
         min_abundance    = 0, 
         min_prevalence   = 0,
         cores            = 1)

Maaslin2(input_data = microbiomes, 
         input_metadata = ddr_plot[,c('participant','a_p22_5_fs_dr')],
         output           = './noGit/maaslin2/DPA_a_p22_5_fs_dr/',
         normalization    = 'NONE', 
         standardize      = 'FALSE',
         transform        = 'AST', 
         analysis_method  = 'LM', 
         random_effects   = 'participant',
         min_abundance    = 0, 
         min_prevalence   = 0,
         cores            = 1)

Maaslin2(input_data = microbiomes, 
         input_metadata = ddr_plot[,c('participant','a_ala_fs_dr')],
         output           = './noGit/maaslin2/ALA_a_ala_fs_dr/',
         normalization    = 'NONE', 
         standardize      = 'FALSE',
         transform        = 'AST', 
         analysis_method  = 'LM', 
         random_effects   = 'participant',
         min_abundance    = 0, 
         min_prevalence   = 0,
         cores            = 1)

Maaslin2(input_data = microbiomes, 
         input_metadata = ddr_plot[,c('participant','a_omega3_noala_fs_dr')],
         output           = './noGit/maaslin2/OMEGA3woala_a_omega3_noala_fs_dr/',
         normalization    = 'NONE', 
         standardize      = 'FALSE',
         transform        = 'AST', 
         analysis_method  = 'LM', 
         random_effects   = 'participant',
         min_abundance    = 0, 
         min_prevalence   = 0,
         cores            = 1)

Maaslin2(input_data = microbiomes, 
         input_metadata = ddr_plot[,c('participant','a_omega3_fs_dr')],
         output           = './noGit/maaslin2/OMEGA3wala_a_a_omega3_fs_dr/',
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
  plot
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

###*****after meeting ****##

scatter_plot_negative('a_omega3_noala_fs_dr','eubacterium_sp_cag_251',"a_omega3_noala_fs_dr\n(g/d)","eubacterium_sp_cag_251\n(relative abundance")
summary(lm(a_omega3_noala_fs_dr ~eubacterium_sp_cag_251, data = ddr_plot))
# Multiple R-squared:  0.002829,	Adjusted R-squared:  0.001728 
# F-statistic:  2.57 on 1 and 906 DF,  p-value: 0.1092
summary(lm(a_omega3_noala_fs_dr ~ eubacterium_sp_cag_251 + agemlvs +abx +calor_fs_dr +(1|participant), data = ddr_plot))
# Multiple R-squared:  0.01016,	Adjusted R-squared:  0.005762 
# F-statistic: 2.311 on 4 and 901 DF,  p-value: 0.0561


scatter_plot_negative('a_omega3_noala_fs_dr','odoribacter_splanchnicus',"a_omega3_noala_fs_dr\n(g/d)","odoribacter_splanchnicus\n(relative abundance")
summary(lm(a_omega3_noala_fs_dr ~ odoribacter_splanchnicus + agemlvs +abx +calor_fs_dr +(1 | participant), data = ddr_plot))
# Multiple R-squared:  0.01932,	Adjusted R-squared:  0.01496 
# F-statistic: 4.437 on 4 and 901 DF,  p-value: 0.00148
summary(lm(a_omega3_noala_fs_dr ~ odoribacter_splanchnicus + agemlvs +abx +calor_fs_dr , data = ddr_plot))
# Multiple R-squared:  0.01932,	Adjusted R-squared:  0.01496 
# F-statistic: 4.437 on 4 and 901 DF,  p-value: 0.00148



scatter_plot_positive('a_omega3_noala_fs_dr','enterorhabdus_caecimuris',"a_omega3_noala_fs_dr\n(g/d)","enterorhabdus_caecimuris\n(relative abundance")
summary(lm(a_omega3_noala_fs_dr ~ enterorhabdus_caecimuris + agemlvs +abx +calor_fs_dr +(1 | participant), data = ddr_plot))
# Multiple R-squared:  0.03004,	Adjusted R-squared:  0.02573 
# F-statistic: 6.976 on 4 and 901 DF,  p-value: 0.00001566
summary(lm(a_omega3_noala_fs_dr ~ enterorhabdus_caecimuris + agemlvs +abx +calor_fs_dr , data = ddr_plot))



#####EPA

scatter_plot_negative('a_f205_fs_dr','eubacterium_sp_cag_251',"a_f205_fs_dr\n(g/d)","eubacterium_sp_cag_251\n(relative abundance")
summary(lm(a_f205_fs_dr ~ eubacterium_sp_cag_251 + agemlvs +abx +calor_fs_dr +(1|participant), data = ddr_plot))
# Multiple R-squared:  0.01037,	Adjusted R-squared:  0.005976 
# F-statistic:  2.36 on 4 and 901 DF,  p-value: 0.0518


scatter_plot_negative('a_f205_fs_dr','odoribacter_splanchnicus',"a_f205_fs_dr\n(g/d)","odoribacter_splanchnicus\n(relative abundance")
summary(lm(a_f205_fs_dr ~ odoribacter_splanchnicus + agemlvs +abx +calor_fs_dr +(1 | participant), data = ddr_plot))
# Multiple R-squared:  0.0184,	Adjusted R-squared:  0.01404 
# F-statistic: 4.221 on 4 and 901 DF,  p-value: 0.002164



scatter_plot_positive('a_f205_fs_dr','enterorhabdus_caecimuris',"a_f205_fs_dr\n(g/d)","enterorhabdus_caecimuris\n(relative abundance")
summary(lm(a_f205_fs_dr ~ enterorhabdus_caecimuris + agemlvs +abx +calor_fs_dr +(1 | participant), data = ddr_plot))
# Multiple R-squared:  0.0388,	Adjusted R-squared:  0.03454 
# F-statistic: 9.094 on 4 and 901 DF,  p-value: 0.0000003337


#####DHA

scatter_plot_negative('a_f226_fs_dr','eubacterium_sp_cag_251',"a_f226_fs_dr\n(g/d)","eubacterium_sp_cag_251\n(relative abundance")
summary(lm(a_f226_fs_dr ~ eubacterium_sp_cag_251 + agemlvs +abx +calor_fs_dr +(1|participant), data = ddr_plot))

scatter_plot_negative('a_f226_fs_dr','odoribacter_splanchnicus',"a_f226_fs_dr\n(g/d)","odoribacter_splanchnicus\n(relative abundance")
summary(lm(a_f226_fs_dr ~ odoribacter_splanchnicus + agemlvs +abx +calor_fs_dr +(1 | participant), data = ddr_plot))

scatter_plot_positive('a_f226_fs_dr','enterorhabdus_caecimuris',"a_f226_fs_dr\n(g/d)","enterorhabdus_caecimuris\n(relative abundance")
summary(lm(a_f226_fs_dr ~ enterorhabdus_caecimuris + agemlvs +abx +calor_fs_dr +(1 | participant), data = ddr_plot))


#####DPA

scatter_plot_negative('a_p22_5_fs_dr','eubacterium_sp_cag_251',"a_p22_5_fs_dr\n(g/d)","eubacterium_sp_cag_251\n(relative abundance")
summary(lm(a_p22_5_fs_dr ~ eubacterium_sp_cag_251 + agemlvs +abx +calor_fs_dr +(1|participant), data = ddr_plot))

scatter_plot_negative('a_p22_5_fs_dr','odoribacter_splanchnicus',"a_p22_5_fs_dr\n(g/d)","odoribacter_splanchnicus\n(relative abundance")
summary(lm(a_p22_5_fs_dr ~ odoribacter_splanchnicus + agemlvs +abx +calor_fs_dr +(1 | participant), data = ddr_plot))

scatter_plot_positive('a_p22_5_fs_dr','enterorhabdus_caecimuris',"a_p22_5_fs_dr\n(g/d)","enterorhabdus_caecimuris\n(relative abundance")
summary(lm(a_p22_5_fs_dr ~ enterorhabdus_caecimuris + agemlvs +abx +calor_fs_dr +(1 | participant), data = ddr_plot))

#####ALA

scatter_plot_positive('a_ala_fs_dr','eubacterium_sp_cag_251',"a_ala_fs_dr\n(g/d)","eubacterium_sp_cag_251\n(relative abundance")
summary(lm(a_ala_fs_dr ~ eubacterium_sp_cag_251 + agemlvs +abx +calor_fs_dr +(1|participant), data = ddr_plot))

scatter_plot_positive('a_ala_fs_dr','odoribacter_splanchnicus',"a_ala_fs_dr\n(g/d)","odoribacter_splanchnicus\n(relative abundance")
summary(lm(a_ala_fs_dr ~ odoribacter_splanchnicus + agemlvs +abx +calor_fs_dr +(1 | participant), data = ddr_plot))

scatter_plot_positive('a_ala_fs_dr','enterorhabdus_caecimuris',"a_ala_fs_dr\n(g/d)","enterorhabdus_caecimuris\n(relative abundance")
summary(lm(a_ala_fs_dr ~ enterorhabdus_caecimuris + agemlvs +abx +calor_fs_dr +(1 | participant), data = ddr_plot))




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
  plot
}# end of function

scatter_plot_negative(species_metadata$a_p22_5_fs_dr, species_metadata$s__Eubacterium_sp_CAG_251, 
                      "a_p22_5_fs_dr", "s__Eubacterium_sp_CAG_251")

species_metadata$a_p22_5_fs_dr
ggplot(data=species_metadata, aes(x=s__Eubacterium_sp_CAG_251, y= a_p22_5_fs_dr)) +
  geom_point(aes(),fill = '#C92D39', color = 'black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
  stat_smooth(method = "glm", color ='black')+ 
  guides(alpha='none')+labs("")+
 # xlab("a_p22_5_fs_dr") +  ylab(ylab)  +
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

species_metadata$s__Eubacterium
species_metadata <-
  species_metadata%>%
  mutate(bi_s__Eubacterium_sp_CAG_251 = ifelse(s__Eubacterium_sp_CAG_251==0, 0, 1))

lmer(data=species_metadata, s__Eubacterium_sp_CAG_251  )

ggplot(data=species_metadata, aes(y=s__Odoribacter_splanchnicus, x= a_omega3_noala_fs_dr)) +
  geom_point(aes(),fill = '#C92D39', color = 'black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
  stat_smooth(method = "glm", color ='black')+ 
  guides(alpha='none')+labs("")+
  # xlab("a_p22_5_fs_dr") +  ylab(ylab)  +
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

species_metadata$s__Odoribacter_splanchnicus
ggplot(data=species_metadata, aes(x= logcrp, y=s__Odoribacter_splanchnicus)) +
  geom_point(aes(),fill = '#C92D39', color = 'black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
  stat_smooth(method = "glm", color ='black')+ 
  guides(alpha='none')+labs("")+
  # xlab("a_p22_5_fs_dr") +  ylab(ylab)  +
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

species_metadata$a_omega3_noala_fs_dr

ggplot(data=species_metadata, aes(x= omega3_noala10v, y=s__Paraprevotella_clara)) +
  geom_point(aes(),fill = '#C92D39', color = 'black', alpha = .5, shape = 21, size = 3, stroke = 1) + 
  stat_smooth(method = "glm", color ='black')+ 
  guides(alpha='none')+labs("")+
  # xlab("a_p22_5_fs_dr") +  ylab(ylab)  +
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




ggplot(data=species_metadata, aes(x=a_p22_5_fs_dr, y=s__Eubacterium_sp_CAG_251))+
  geom_point()

species_metadata <-
  species_metadata%>%
  mutate(bi_s__Odoribacter_splanchnicus = ifelse(s__Odoribacter_splanchnicus==0, 0, 1))

summary(lmer(data = species_metadata, logcrp ~ a_omega3_noala_fs_dr + s__Odoribacter_splanchnicus 
                                        +a_omega3_noala_fs_dr:s__Odoribacter_splanchnicus
                                        +agemlvs +abx +calor_fs_dr +(1|participant)))

summary(lmer(data = species_metadata, logcrp ~ a_omega3_noala_fs_dr + bi_s__Eubacterium_sp_CAG_251 
             +a_omega3_noala_fs_dr:bi_s__Eubacterium_sp_CAG_251
             +agemlvs +abx +calor_fs_dr +(1|participant)))
confint(lmer(data = species_metadata, logcrp ~ a_omega3_noala_fs_dr + bi_s__Eubacterium_sp_CAG_251 
             +a_omega3_noala_fs_dr:bi_s__Eubacterium_sp_CAG_251
             +agemlvs +abx +calor_fs_dr +(1|participant)))


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


scatter_plot_negative('logcrp','eubacterium_eligens',"Short-term dietary fiber\n(g/d)","Eubacterium eligens\n(relative abundance")












# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
####################### Fig 3A. MaAsLin2 ################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################# read in data ###################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

pfa_file <- c('ala_pfa', 'epa_pfa', 'dpa_pfa', 'dha_pfa')

Maaslin2(input_data = microbiomes, 
         input_metadata = cbind(adjust, ala_pfa = ddr_plot[,c('ala_pfa')]),
         output           = './noGit/maaslin2/pfa_ala/',
         normalization    = 'NONE', 
         standardize      = 'FALSE',
         transform        = 'AST', 
         analysis_method  = 'LM', 
         random_effects   = 'participant',
         min_abundance    = 0, 
         min_prevalence   = 0,
         cores            = 1)

Maaslin2(input_data = microbiomes, 
         input_metadata = cbind(adjust, epa_pfa = ddr_plot[,c('epa_pfa')]),
         output           = './noGit/maaslin2/pfa_epa/',
         normalization    = 'NONE', 
         standardize      = 'FALSE',
         transform        = 'AST', 
         analysis_method  = 'LM', 
         random_effects   = 'participant',
         min_abundance    = 0, 
         min_prevalence   = 0,
         cores            = 1)

Maaslin2(input_data = microbiomes, 
         input_metadata = cbind(adjust, dha_pfa = ddr_plot[,c('dha_pfa')]),
         output           = './noGit/maaslin2/pfa_dha/',
         normalization    = 'NONE', 
         standardize      = 'FALSE',
         transform        = 'AST', 
         analysis_method  = 'LM', 
         random_effects   = 'participant',
         min_abundance    = 0, 
         min_prevalence   = 0,
         cores            = 1)

Maaslin2(input_data = microbiomes, 
         input_metadata = cbind(adjust, dpa_pfa = ddr_plot[,c('dpa_pfa')]),
         output           = './noGit/maaslin2/pfa_dpa/',
         normalization    = 'NONE', 
         standardize      = 'FALSE',
         transform        = 'AST', 
         analysis_method  = 'LM', 
         random_effects   = 'participant',
         min_abundance    = 0, 
         min_prevalence   = 0,
         cores            = 1)




pfa_file <- c('ala_pfa', 'epa_pfa', 'dpa_pfa', 'dha_pfa')


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################## @@ Figure 3B/C. scatter plots ###################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

######### @ need modification from fiber to omega3# ########

# Function of positive associations

#####EPA

scatter_plot_positive('epa_pfa','firmicutes_bacterium_cag_94',"epa_pfa\n(g/d)","firmicutes_bacterium_cag_94\n(relative abundance")
summary(lm(epa_pfa ~ firmicutes_bacterium_cag_94 + agemlvs +abx +calor_fs_dr +(1|participant), data = ddr_plot))
# Multiple R-squared:  0.04139,	Adjusted R-squared:  0.03703 
# F-statistic: 9.509 on 4 and 881 DF,  p-value: 0.0000001578


scatter_plot_negative('epa_pfa','fusicatenibacter_saccharivorans',"epa_pfa\n(g/d)","fusicatenibacter_saccharivorans\n(relative abundance")
summary(lm(epa_pfa ~ fusicatenibacter_saccharivorans + agemlvs +abx +calor_fs_dr +(1 | participant), data = ddr_plot))
# Multiple R-squared:  0.01982,	Adjusted R-squared:  0.01537 
# F-statistic: 4.454 on 4 and 881 DF,  p-value: 0.001438



#####EPA

scatter_plot_negative('epa_pfa','eubacterium_sp_cag_251',"epa_pfa\n(g/d)","eubacterium_sp_cag_251\n(relative abundance")
summary(lm(epa_pfa ~ eubacterium_sp_cag_251 + agemlvs +abx +calor_fs_dr +(1|participant), data = ddr_plot))
# Multiple R-squared:  0.009069,	Adjusted R-squared:  0.00457 
# F-statistic: 2.016 on 4 and 881 DF,  p-value: 0.09029


scatter_plot_negative('epa_pfa','odoribacter_splanchnicus',"epa_pfa\n(g/d)","odoribacter_splanchnicus\n(relative abundance")
summary(lm(epa_pfa ~ odoribacter_splanchnicus + agemlvs +abx +calor_fs_dr +(1 | participant), data = ddr_plot))
# # Multiple R-squared:  0.01594,	Adjusted R-squared:  0.01148 
# F-statistic: 3.569 on 4 and 881 DF,  p-value: 0.006755



scatter_plot_positive('epa_pfa','enterorhabdus_caecimuris',"epa_pfa\n(g/d)","enterorhabdus_caecimuris\n(relative abundance")
summary(lm(epa_pfa ~ enterorhabdus_caecimuris + agemlvs +abx +calor_fs_dr +(1 | participant), data = ddr_plot))
# Multiple R-squared:  0.0535,	Adjusted R-squared:  0.0492 
# F-statistic: 12.45 on 4 and 881 DF,  p-value: 0.0000000007445






#####ALA

scatter_plot_positive('a_ala_fs_dr','eubacterium_sp_cag_251',"a_ala_fs_dr\n(g/d)","eubacterium_sp_cag_251\n(relative abundance")
summary(lm(a_ala_fs_dr ~ eubacterium_sp_cag_251 + agemlvs +abx +calor_fs_dr +(1|participant), data = ddr_plot))

scatter_plot_positive('a_ala_fs_dr','odoribacter_splanchnicus',"a_ala_fs_dr\n(g/d)","odoribacter_splanchnicus\n(relative abundance")
summary(lm(a_ala_fs_dr ~ odoribacter_splanchnicus + agemlvs +abx +calor_fs_dr +(1 | participant), data = ddr_plot))

scatter_plot_positive('a_ala_fs_dr','enterorhabdus_caecimuris',"a_ala_fs_dr\n(g/d)","enterorhabdus_caecimuris\n(relative abundance")
summary(lm(a_ala_fs_dr ~ enterorhabdus_caecimuris + agemlvs +abx +calor_fs_dr +(1 | participant), data = ddr_plot))




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#######################  intake and plasma concentration ################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
################# read in data ###################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#####EPA

scatter_plot_positive('a_f205_fs_dr','epa_pfa',"a_f205_fs_dr\n(g/d)","epa_pfa\n(relative abundance")
summary(lm(a_f205_fs_dr ~ epa_pfa + agemlvs +abx +calor_fs_dr +(1|participant), data = z_mlvs_exposure_plasma))
# Multiple R-squared:  0.04139,	Adjusted R-squared:  0.03703 
# F-statistic: 9.509 on 4 and 881 DF,  p-value: 0.0000001578


scatter_plot_negative('epa_pfa','fusicatenibacter_saccharivorans',"epa_pfa\n(g/d)","fusicatenibacter_saccharivorans\n(relative abundance")
summary(lm(epa_pfa ~ fusicatenibacter_saccharivorans + agemlvs +abx +calor_fs_dr +(1 | participant), data = ddr_plot))
# Multiple R-squared:  0.01982,	Adjusted R-squared:  0.01537 
# F-statistic: 4.454 on 4 and 881 DF,  p-value: 0.001438



#####EPA

scatter_plot_negative('epa_pfa','eubacterium_sp_cag_251',"epa_pfa\n(g/d)","eubacterium_sp_cag_251\n(relative abundance")
summary(lm(epa_pfa ~ eubacterium_sp_cag_251 + agemlvs +abx +calor_fs_dr +(1|participant), data = ddr_plot))
# Multiple R-squared:  0.009069,	Adjusted R-squared:  0.00457 
# F-statistic: 2.016 on 4 and 881 DF,  p-value: 0.09029


scatter_plot_negative('epa_pfa','odoribacter_splanchnicus',"epa_pfa\n(g/d)","odoribacter_splanchnicus\n(relative abundance")
summary(lm(epa_pfa ~ odoribacter_splanchnicus + agemlvs +abx +calor_fs_dr +(1 | participant), data = ddr_plot))
# # Multiple R-squared:  0.01594,	Adjusted R-squared:  0.01148 
# F-statistic: 3.569 on 4 and 881 DF,  p-value: 0.006755



scatter_plot_positive('epa_pfa','enterorhabdus_caecimuris',"epa_pfa\n(g/d)","enterorhabdus_caecimuris\n(relative abundance")
summary(lm(epa_pfa ~ enterorhabdus_caecimuris + agemlvs +abx +calor_fs_dr +(1 | participant), data = ddr_plot))
# Multiple R-squared:  0.0535,	Adjusted R-squared:  0.0492 
# F-statistic: 12.45 on 4 and 881 DF,  p-value: 0.0000000007445






#####ALA

scatter_plot_positive('a_ala_fs_dr','eubacterium_sp_cag_251',"a_ala_fs_dr\n(g/d)","eubacterium_sp_cag_251\n(relative abundance")
summary(lm(a_ala_fs_dr ~ eubacterium_sp_cag_251 + agemlvs +abx +calor_fs_dr +(1|participant), data = ddr_plot))

scatter_plot_positive('a_ala_fs_dr','odoribacter_splanchnicus',"a_ala_fs_dr\n(g/d)","odoribacter_splanchnicus\n(relative abundance")
summary(lm(a_ala_fs_dr ~ odoribacter_splanchnicus + agemlvs +abx +calor_fs_dr +(1 | participant), data = ddr_plot))

scatter_plot_positive('a_ala_fs_dr','enterorhabdus_caecimuris',"a_ala_fs_dr\n(g/d)","enterorhabdus_caecimuris\n(relative abundance")
summary(lm(a_ala_fs_dr ~ enterorhabdus_caecimuris + agemlvs +abx +calor_fs_dr +(1 | participant), data = ddr_plot))






##### draw Fig 3 #####
maaslin <- read.csv("noGit/maaslin.csv")
signifbugs <- maaslin$feature
signifbugs <- signifbugs[!duplicated(signifbugs)]
signifbugs

signifmeta <- maaslin$metadata
signifmeta <- signifmeta[!duplicated(signifmeta)]



ggplot(data = maaslin, aes(x=metadata, y=feature))+
  geom_tile(aes(fill = coef ), colour = "black")+
  scale_fill_gradientn(colours = c( "#1207A3", "#1207A3", "#FFFEFE","#BB0103"))+# "#BB0103" "#8F88D2", "#FFFEFE","#DB7A7B"
  geom_text(aes(label = ifelse(pval<0.1, "*", NA)) )+
  theme_classic()+
  theme(axis.text.x=element_text(angle=45,hjust=1))

maaslin2 <- read.csv("noGit/maaslin2.csv")
summary(maaslin2$metadata)
maaslin2$metadata <- factor(maaslin2$metadata, 
                               levels=c("ala10v",  "epa10v", "dha10v", "omega310v", 
                                        "ala_avg", "epa_avg",  "dha_avg","dpa_avg",  "omega3_avg", "omega3_noala_avg",  "omega6_avg", "trans_avg",
                                        "a_ala_fs_dr", "a_f205_fs_dr",  "a_f226_fs_dr",  "a_p22_5_fs_dr",  "a_omega3_fs_dr", "a_omega3_noala_fs_dr",  "a_trn07_fo_dr",    
                                        "logcrp",  "loghdl", "logtg" ) )
metadd <- c("ala10v",  "epa10v", "dha10v", "omega310v", 
            "ala_avg", "epa_avg",  "dha_avg","dpa_avg",  "omega3_avg", "omega3_noala_avg",  "omega6_avg", "trans_avg",
            "a_ala_fs_dr", "a_f205_fs_dr",  "a_f226_fs_dr",  "a_p22_5_fs_dr",  "a_omega3_fs_dr", "a_omega3_noala_fs_dr",  "a_trn07_fo_dr",    
            "logcrp",  "loghdl", "logtg" )                    
 

ggplot(data = maaslin2, aes(x=metadata, y=feature))+
  geom_tile(aes(fill = ifelse(coef>0.1, 0.1, 
                              ifelse(coef< -0.1, -0.1, coef))), colour = "black")+
  scale_fill_gradientn(name = "coef, *qval<0.10", colours = c( "#1207A3", "#FFFEFE","#BB0103"),
                       )+# "#BB0103" "#8F88D2", "#FFFEFE","#DB7A7B"
  geom_text(aes(label = ifelse(qval<0.1, "*", NA)) )+
 # theme_classic()+
  theme(axis.text.x=element_text(angle=45,hjust=1)
        )






mamamama <- rbind(read_tsv("noGit/maaslin2/omega_output_ala_ffqcum_species_all/all_results.tsv") ,
                   read_tsv("noGit/maaslin2/omega_output_epa_ffqcum_species_all/all_results.tsv"),
                   read_tsv("noGit/maaslin2/omega_output_dha_ffqcum_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_dpa_ffqcum_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_trans_ffqcum_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_omega6_ffqcum_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_omega3_ffqcum_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_omega3_noala_ffqcum_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_ala_avg_fs_ffq_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_epa_avg_fs_ffq_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_dha_avg_fs_ffq_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_dpa_avg_fs_ffq_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_trans_avg_fs_ffq_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_omega6_avg_fs_ffq_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_omega3_avg_fs_ffq_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_omega3_noala_avg_fs_ffq_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_a_omega3_fs_ddr_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_a_trans_fs_ddr_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_a_ala_fs_ddr_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_a_epa_fs_ddr_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_a_dha_fs_ddr_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_a_dpa_fs_ddr_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_a_omega3_noala_fs_ddr_species_all/all_results.tsv"),
                  
                  read_tsv("noGit/maaslin2/omega_output_biomarker_logcrp_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_biomarker_loghdl_species_all/all_results.tsv"),
                  read_tsv("noGit/maaslin2/omega_output_biomarker_logtg_species_all/all_results.tsv"))

mam <- mamamama %>%
  filter(metadata == metadd) #%>%
  filter(feature == signifbugs)

ggplot(data = mamamama, aes(x=metadata, y=feature))+
  geom_tile(aes(fill = coef ), colour = "black")+
  scale_fill_gradientn(colours = c( "#1207A3", "#1207A3", "#FFFEFE","#BB0103"))+# "#BB0103" "#8F88D2", "#FFFEFE","#DB7A7B"
  geom_text(aes(label = ifelse(pval<0.1, "*", NA)) )+
  theme_classic()+
  theme(axis.text.x=element_text(angle=45,hjust=1))








  summary(species_metadata$a_f205_fs_dr)
  ggplot(species_metadata, aes(x=a_f205_fs_dr))+ #a_omega3_fs_dr_w1
    geom_histogram(color="black", fill="grey")+#, binwidth=0.5)+
  labs(x="a_f205_fs_dr(EPA) (g/day)", y = "Count")+
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

  
  species_metadata$s__Eubacterium_sp_CAG_251
  
  library(lme4)
  
  hithere <- species_metadata[c('a_p22_5_fs_dr', 's__Eubacterium_sp_CAG_251', 
                                'participant', 'agemlvs', 'calor_fs_dr', 'abx')]
  hithere <- hithere%>%
    mutate(bi_s__Eubacterium_sp_CAG_251 = ifelse(s__Eubacterium_sp_CAG_251==0, 0, 1))
  
  summary(lmer(data = hithere, a_p22_5_fs_dr ~bi_s__Eubacterium_sp_CAG_251 +agemlvs+ calor_fs_dr+ abx+(1|participant)))
library(lmerTest)  
  
  