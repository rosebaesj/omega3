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
                    'biomarker_logcrp', 'biomarker_logcrp_complete', 'a_omega3_fs_ddr_logcrp', 'biomarker_loghdl', 'biomarker_loghdl_complete', 'a_omega3_fs_ddr_loghdl', 'biomarker_logtg', 'biomarker_logtg_complete', 'a_omega3_fs_ddr_logtg', 'a_trans_fs_ddr_logcrp', 'a_trans_fs_ddr_loghdl', 'a_trans_fs_ddr_logtg',
                    'a_ala_fs_ddr_logcrp', 'a_ala_fs_ddr_loghdl', 'a_ala_fs_ddr_logtg',
                    'a_epa_fs_ddr_logcrp', 'a_epa_fs_ddr_loghdl', 'a_epa_fs_ddr_logtg',
                    'a_dha_fs_ddr_logcrp', 'a_dha_fs_ddr_loghdl', 'a_dha_fs_ddr_logtg',
                    'a_dpa_fs_ddr_logcrp', 'a_dpa_fs_ddr_loghdl', 'a_dpa_fs_ddr_logtg')

pfa_file <- c('ala_pfa', 'epa_pfa', 'dpa_pfa', 'dha_pfa')

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


