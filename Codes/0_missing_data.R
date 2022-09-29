##############################################################################################################################################
#+ 1) Purposes: conduct correlation analysis among the main exposure and outcome variables, 
#+              Main exposure: omega3 - long/short term, ala epa dha dpa
#+ 2) Study design:  Cross-sectional design
#+ 3) Endpoints: primary: 
# 4) Exposures: primary: cumulative average of long-term physical activity level
# 5) Covariates: age, diet quality, total energy intake, smoking status, antibiotic use, probiotic use, stool type
# 6) Follow-up: 2012 (Men's Lifestyle Validation Study)  
#############################################################################################################################################




###### DATA PREPARATION ########

####### 1) Microbiome Taxonomic Table ######

#### TAXONOMY DATA ####
# read in taxonomy data
# only keep species-level features

all_species <- read_tsv("data/metaphlan_taxonomic_profiles.tsv")
names(all_species)[names(all_species) == '# taxonomy'] <- 'Sample'

all_species <-all_species %>%
  separate(Sample, c("kingdom","phylum","class" ,"order","family","genus" ,"species" ), 
           sep = '\\|', remove = TRUE)

all_species1 <- subset(all_species,!is.na(species))
all_species1 <- as.data.frame(all_species1) #weird formatting change to maintain rows
rownames(all_species1)<-all_species1$species

###### TAX and OTU table
TAX <- as.matrix(all_species1[,1:7])
OTU <- all_species1[,8:ncol(all_species1)]


all_species1<-all_species1[,-c(1:7)]
all_species1 <- as.data.frame(all_species1)
dim(all_species1)
# sunjeong: [1] 543 929

all_sample <- colnames(all_species1)
#View(all_sample)

# remove guy with colectomy 
all_species1 <- all_species1[ , !colnames(all_species1) %in% c(grep("15074981", colnames(all_species1), value = T))]

# check NA
sum(colSums(is.na(all_species1)))
#0 -> no 



# generate list of species after abundance/prevalence filtering
# since this file is scaled to 100, i am keeping bugs with >10% prevalence and 0.01 abundance (i.e. .0001*100 when relab is scaled 0-1 and not 0-100)
dim(all_species1[apply(all_species1, 1, function(bug) sum(bug >= 0.01) >= 0.1*ncol(all_species1)), ])
# [1] 139 925 = the 139 species after QC included in initial starr manuscripts
# sunjeong: [1] 156 925
all_species_filt <- all_species1[ apply(all_species1, 1, function(bug) sum(bug >= 0.01) >= 0.1*ncol(all_species1)), ]


# important! bug tables from channing need to be converted to 0-1 (not 0-11), which will wreck our downstream stats
all_species_filt   <- all_species_filt  / 100 

all_species_filt <- t(all_species_filt)
all_species_filt <- as.data.frame(all_species_filt)

all_species_filt$ID1 <- substr(rownames(all_species_filt), start = 1, stop = 6)
all_species_filt$stn <- substr(rownames(all_species_filt), start = 10, stop = 13)
all_species_filt$stn <- str_replace_all(all_species_filt$stn,  c("SF05" = "s1","SF06" = "s2","SF07" = "s3","SF08" = "s4"))
rownames(all_species_filt) <- paste0(all_species_filt$ID1,"_", all_species_filt$stn)

species_filt <- all_species_filt[,!colnames(all_species_filt) %in% c('ID1','stn')]

##### save the metadata of exposure and plasma outcome ####

write.table(species_filt, file = "./data/species_filt.csv", sep = ",", col.names = TRUE, row.names = TRUE)








### 2) Metadata of Dietary #########

#### + Import Metadata by individual #####
meta_id <- read.table("data/mlvs_exposure_for_jorick.txt", header=TRUE, sep='\t', check.names=TRUE, quote ="")
idkey <- read.csv('data/idkey.csv', header=TRUE, sep = ",", stringsAsFactors = FALSE)
rownames(idkey) <- idkey$id
meta_id<- left_join(idkey, meta_id, by="id")

# missing1 <- is.na(meta_id[,c('logcrp_plasma1', 'a_omega3_fs_dr_w1', 'a_f205_fs_dr_w1', 'a_f226_fs_dr_w1', 'a_p22_5_fs_dr_w1',
#                             'omega3_noala_ffq1', 'epa_ffq1', 'dpa_ffq1', 'dha_ffq1',
#                             'omega3_noala10v', 'epa10v', 'dpa10v', 'dha10v'
#                             ) ])
# # colSums(missing1)
# sum(rowSums(missing1)!=0)
# missing2 <- is.na(meta_id[,c('logcrp_plasma2', 'a_omega3_fs_dr_w2', 'a_f205_fs_dr_w2', 'a_f226_fs_dr_w2', 'a_p22_5_fs_dr_w2',
#                              'omega3_noala_ffq2', 'epa_ffq2', 'dpa_ffq2', 'dha_ffq2',
#                              'omega3_noala10v', 'epa10v', 'dpa10v', 'dha10v') ])
# colSums(missing2)
# sum(rowSums(missing2)!=0)
# 
# missing <- cbind(missing1, missing2)
# colSums(missing)
# sum(rowSums(missing)!=0)
# 
# colSums(is.na(meta_id[,c('logcrp_plasma1', 'a_omega3_fs_dr_w1', 'a_f205_fs_dr_w1', 'a_f226_fs_dr_w1', 'a_p22_5_fs_dr_w1') ]))
# sum(is.na(meta_id$logcrp_plasma1 )|is.na(meta_id$a_omega3_fs_dr_w1)|is.na(meta_id$a_f205_fs_dr_w1)|is.na(meta_id$a_f226_fs_dr_w1)|is.na(meta_id$a_p22_5_fs_dr_w1))
# 
# colSums(is.na(meta_id[,c('logcrp_plasma2', 'a_omega3_fs_dr_w2', 'a_f205_fs_dr_w2', 'a_f226_fs_dr_w2', 'a_p22_5_fs_dr_w2') ]))
# sum(is.na(meta_id$logcrp_plasma2 )|is.na(meta_id$a_omega3_fs_dr_w2)|is.na(meta_id$a_f205_fs_dr_w2)|is.na(meta_id$a_f226_fs_dr_w2)|is.na(meta_id$a_p22_5_fs_dr_w2))
# 
# colSums(is.na(meta_id[,c('logcrp_plasma1', 'omega3_noala_ffq1', 'epa_ffq1', 'dpa_ffq1', 'dha_ffq1') ]))
# sum(is.na(meta_id$logcrp_plasma1 )|is.na(meta_id$omega3_noala_ffq1)|is.na(meta_id$epa_ffq1)|is.na(meta_id$dpa_ffq1)|is.na(meta_id$dha_ffq1))
# 
# colSums(is.na(meta_id[,c('logcrp_plasma1', 'a_omega3_fs_dr_w2', 'epa_ffq1', 'dpa_ffq1', 'dha_ffq1') ]))
# sum(is.na(meta_id$logcrp_plasma1 )|is.na(meta_id$omega3_noala10v)|is.na(meta_id$epa_ffq1)|is.na(meta_id$dpa_ffq1)|is.na(meta_id$dha_ffq1))






#### + set Phase #####
#+ ph0 : long term
#+ ph1 : first two stools
#+ ph2 : second two stools

#### all together ####
colnames(meta_id) <- gsub("w1avg","w1", colnames(meta_id)) 
colnames(meta_id) <- gsub("w2avg","w2", colnames(meta_id))

phase1 <- meta_id %>% select(ID1, contains('_w1'), contains('ffq1'), contains('plasma1'), contains('qu1'))
phase2 <- meta_id %>% select(ID1, contains('_w2'), contains('ffq2'), contains('plasma2'), contains('qu2'))

colnames(phase1) <- str_replace_all(colnames(phase1),  c("_w1" = "_ddr", "ffq1" = "ffq", "plasma1" = "plasma", "qu1" = "qu"))
colnames(phase2) <- str_replace_all(colnames(phase2),  c("_w2" = "_ddr", "ffq2" = "ffq", "plasma2" = "plasma", "qu2" = "qu"))

##### missing data ####
dim(phase1)==dim(phase2)
 
 for (i in 1:ncol(phase1)){
   for (j in 1:nrow(phase1)){
     if (is.na(phase1[j,i])){phase1[j,i]<- phase2[j,i]}
   }
 }
 for (i in 1:ncol(phase2)){
   for (j in 1:nrow(phase2)){
     if (is.na(phase2[j,i])){phase2[j,i]<- phase1[j,i]}
   }
 }

phase1$phase <- c("ph1")
phase2$phase <- c("ph2")
meta_id_nophase <- meta_id %>% select(everything(), 
                                      -contains('_w1'), -contains('ffq1'), -contains('plasma1'), -contains('qu1'),
                                      -contains('_w2'), -contains('ffq2'), -contains('plasma2'), -contains('qu2'))


# #### + a) DDR #####
# colnames(meta_id) <- gsub("w1avg","w1", colnames(meta_id)) 
# colnames(meta_id) <- gsub("w2avg","w2", colnames(meta_id))
# 
# ddr1 <- meta_id %>% select(ID1, contains('_w1'))
# ddr2 <- meta_id %>% select(ID1, contains('_w2'))
# 
# colnames(ddr1) <- str_replace_all(colnames(ddr1),  c("_w1" = "_ddr"))
# colnames(ddr2) <- str_replace_all(colnames(ddr2),  c("_w2" = "_ddr"))
# 
# 
# dim(ddr1)==dim(ddr2)
# 
# for (i in 1:ncol(ddr1)){
#   for (j in 1:nrow(ddr1)){
#     if (is.na(ddr1[j,i])){ddr1[j,i]<- ddr2[j,i]}
#   }
# }
# for (i in 1:ncol(ddr2)){
#   for (j in 1:nrow(ddr2)){
#     if (is.na(ddr2[j,i])){ddr2[j,i]<- ddr1[j,i]}
#   }
# }
# ddr1$phase <- c("ph1")
# ddr2$phase <- c("ph2")
# meta_id_noddr <- meta_id %>% select(everything(), -contains('_w1'), -contains('_w2'))



# #### + e) avg is from ffq1 ffq2, qu1 qu2 ####
# avg <- meta_id_noddr_noffq_noplasma_noqu %>% select(ID1, contains('_avg'))
# 
# meta_id_noddr_noffq_noplasma_noqu_noavg <- meta_id_noddr_noffq_noplasma_noqu %>% select(everything(), -contains('avg'))
# 
# #### + f) 10v is cummulative ffq ####
# cum10v <- meta_id_noddr_noffq_noplasma_noqu_noavg %>% select(ID1, contains('10v'))
# 
# meta_id_noddr_noffq_noplasma_noqu_noavg_no10v <- meta_id_noddr_noffq_noplasma_noqu_noavg %>% select(everything(), -contains('10v'))
# 
# basic <- meta_id_noddr_noffq_noplasma_noqu_noavg_no10v


#### LETS MERGE METADATA with sample num ####

all_species_filt <- all_species_filt %>%
  mutate(phase = ifelse(stn=="s1"|stn=="s2","ph1",
                        ifelse(stn=="s3"|stn=="s4", "ph2", NA))) %>%
  mutate(ID1_stn = rownames(all_species_filt))

all_species_filt$ID1 <- as.integer(all_species_filt$ID1)
samplenum <- all_species_filt %>% select(ID1_stn, ID1) 
samplenum_ph1 <- all_species_filt %>% 
  filter(phase=="ph1")%>%
  select(ID1_stn, ID1)
samplenum_ph2 <- all_species_filt %>% 
  filter(phase=="ph2")%>%
  select(ID1_stn, ID1)


#colnames(phase1)==colnames(phase2)
colnames(phase1)<-colnames(phase2)

meta_id_nophase<- dplyr::rename(meta_id_nophase, pha= phase)

samplenum_noph <- left_join(samplenum, meta_id_nophase, by = "ID1")

samplenum_ph1 <- left_join(samplenum_ph1, phase1, by = "ID1")
samplenum_ph2 <- left_join(samplenum_ph2, phase2, by = "ID1")

samplenum_ph <- rbind(samplenum_ph1, samplenum_ph2)

samplenum_ph <- samplenum_ph %>% select(everything(), -'ID1')

metadata_diet <- left_join(samplenum_noph, samplenum_ph, by = "ID1_stn")


#### 3) Bristol score #####

bristol5 <- metadata_diet %>% 
  mutate(ID1_stn = paste(ID1,"_s1")) %>%
  select(ID1_stn, bristol=bristol5, bristolcat=bristolcat5)

bristol6 <- metadata_diet %>% 
  mutate(ID1_stn = paste(ID1,"_s2")) %>%
  select(ID1_stn, bristol=bristol6, bristolcat=bristolcat6)

bristol7 <- metadata_diet %>% 
  mutate(ID1_stn = paste(ID1,"_s3")) %>%
  select(ID1_stn, bristol=bristol7, bristolcat=bristolcat7)

bristol8 <- metadata_diet %>% 
  mutate(ID1_stn = paste(ID1,"_s4")) %>%
  select(ID1_stn, bristol=bristol8, bristolcat=bristolcat8)

##### missing data ####
# for (i in 1:nrow(bristol5)){
#   if (is.na(bristol5[i,2])){bristol5[i,2]<- metadata_diet$bristol_avg[i]}
#   if (is.na(bristol6[i,2])){bristol6[i,2]<- metadata_diet$bristol_avg[i]}
#   if (is.na(bristol7[i,2])){bristol7[i,2]<- metadata_diet$bristol_avg[i]}
#   if (is.na(bristol8[i,2])){bristol8[i,2]<- metadata_diet$bristol_avg[i]}
#   
#   if (is.na(bristol5[i,3])){bristol5[i,3]<- metadata_diet$bristolcat_avg[i]}
#   if (is.na(bristol6[i,3])){bristol6[i,3]<- metadata_diet$bristolcat_avg[i]}
#   if (is.na(bristol7[i,3])){bristol7[i,3]<- metadata_diet$bristolcat_avg[i]}
#   if (is.na(bristol8[i,3])){bristol8[i,3]<- metadata_diet$bristolcat_avg[i]}
# }

bristol <- rbind(bristol5, bristol6, bristol7, bristol8)
metadata_db <- left_join(metadata_diet, bristol, by='ID1_stn')

metadata_db <- metadata_db %>% select(everything(), -'bristol5', -'bristol6', -'bristol7', -'bristol8',
                                            -'bristolcat5', -'bristolcat6', -'bristolcat7', -'bristolcat8')







#### 3) Import Metabolite data by individual#####
pfa <- read.table("data/plasma_fatty_acids_sjbae.txt",
                     header=TRUE, sep='\t', check.names=TRUE, quote ="")
pfa <- pfa %>%
  mutate(ID1 = as.integer(substr(pfa$ID, 1, 6)))%>%
  mutate(omega3_pfa = ala_pfa+ epa_pfa+ dpa_pfa+ dha_pfa,
         omega3_noala_pfa= epa_pfa+ dpa_pfa+ dha_pfa,
         omega6_pfa = P1+P2+P5+ P6+ P8+P9+P12)


#get first 6 number of ID1 which stands for id in idkey file

#### Merge and get individual metadata 
metadata_dbp<- left_join(metadata_db, pfa, by ="ID1")

metadata_dbp$calor_fs_dr_ddr



### 4) metadata of interest ####
basic <-   c('ID1_stn', 'ID1', 'phase',
           'agemlvs', 'bmi10', 'bmi12', 'wt10', 'wt12', 'smoke10', 'smk12', 'ncig12', 'bmic10', 'bmic12', 'bmi12cat', 'smoke12', 
           'calor10n', 'alco10n', 'calor10v', 
           'calor_avg', 'alco_avg', 'abx_avg', 'probx_avg', 'alc_avg', 'bristol', 'bristolcat', 'bristol_avg', 'bristolcat_avg',
           'a_alco_fo_dr_ddr', 'abx_ddr', 'probx_ddr', 'bristol_ddr', 'bristolcat_ddr', 'calor_fs_dr_ddr'  )

# smoke12 is ordinal
exposure_long <- c('ala10v', 'epa10v', 'dha10v', 'trans10v', 'omega610v', 'omega3_noala10a', 'dpa10v', 'omega310v', 'omega3_noala10v')
exposure_avg <- c('ala_avg', 'epa_avg', 'dpa_avg', 'trans_avg', 'omega3_avg', 'omega6_avg', 'omega3_noala_avg', 'fishavg')
exposure_short <- c( 'ala_ddr' = 'a_ala_fs_dr_ddr',
                     'epa_ddr' = 'a_f205_fs_dr_ddr',
                     'dha_ddr' = 'a_f226_fs_dr_ddr',
                     'dpa_ddr' = 'a_p22_5_fs_dr_ddr',
                     'omega3_noala_ddr' = 'a_omega3_noala_fs_dr_ddr',
                     'omega3_ddr' = 'a_omega3_fs_dr_ddr',
                        'fat' = 'a_fat_fs_dr_ddr')

exposure_1yr <- c('ala_ffq', 'epa_ffq', 'dpa_ffq', 'dha_ffq', 'omega3_ffq', 'omega6_ffq', 'omega3_noala_ffq', 'omega_3_epa_ffq',
                  'tuna_ffq', 'dk_fish_ffq', 'oth_fish_ffq','fishffq')
outcome_plasma <- c('folate_plasma', 'tc_plasma', 'hdlc_plasma', 'tg_plasma', 'crp_plasma',  'crp_plasmacat', 'tchdl_plasma', 
            'logtchdl_plasma', 'logcrp_plasma', 'logtc_plasma', 'loghdl_plasma', 'logtg_plasma')
outcome_pfa <- c('ala_pfa', 'epa_pfa', 'dpa_pfa', 'dha_pfa', 
                'AA_pfa' = 'P8', 
                'omega3_pfa','omega3_noala_pfa','omega6_pfa')
outcome_logpfa <- c('logala_pfa', 'logepa_pfa', 'logdpa_pfa', 'logdha_pfa', 
                    'logAA_pfa' = 'logP8', 
                    'logomega3_pfa','logomega3_noala_pfa','logomega6_pfa')

metadata_intrst <- metadata_dbp %>%
  select(basic, exposure_long, exposure_avg, exposure_1yr, exposure_short,
         outcome_plasma, outcome_pfa)



###3) Missing data ####
# median <- metadata_intrst %>% 
#   dplyr::summarise(across(everything(), ~ mean(.x, na.rm=T))) 
# 
# for (i in 1:ncol(metadata_intrst)){
#   for( j in 1:nrow(metadata_intrst)){
#     if(is.na(metadata_intrst[j,i])){metadata_intrst[j,i]<- median[1,i]}
#   }
# }

#metadata_intrst <- tibble::column_to_rownames(metadata_intrst, var = "ID1_stn")


#### Make log ######
for_log <- metadata_intrst
for_log <- for_log %>% 
  dplyr::select(where(is.numeric))

for_log <- for_log[, colSums(for_log, na.rm = T)!=0]

for (i in 1:ncol(for_log)){
  a<-for_log[,i]
  for (j in 1:nrow(for_log)){
    if (for_log[j,i]<=0){for_log[j,i]<- min(a[a > 0])}
  }
}


logmeta <- log(for_log[,-1])
colnames(logmeta) <- paste("log", colnames(for_log[,-1]), sep = "")
logmeta <- cbind("ID1_stn" = metadata_intrst$ID1_stn, "ID1" = metadata_intrst$ID1, logmeta)

totmeta <- cbind(metadata_intrst, logmeta[,-c(1:2)])






## log species too ####
logspecies <- species_filt

logspecies <- logspecies[, colSums(logspecies, na.rm = T)!=0]

for (i in 1:ncol(logspecies)){
  a<-logspecies[,i]
  for (j in 1:nrow(logspecies)){
    if (logspecies[j,i]<=0){logspecies[j,i]<- min(a[a > 0])}
  }
}


logspecies <- log(logspecies)

## binary species too ####
bispecies <- species_filt

bispecies <- bispecies[, colSums(bispecies, na.rm = T)!=0]

for (i in 1:ncol(bispecies)){
  mm <- median(bispecies[,i])
  for (j in 1:nrow(bispecies)){
    bispecies[j,i] <- ifelse(bispecies[j,i]<=mm, 0, 1)
  }
}




###DONE!!!!!!!!!!!!!!!!!!
missing1 <- is.na(metadata_intrst[,c('logcrp_plasma', 'omega3_noala_ddr', 'epa_ddr', 'dpa_ddr', 'dha_ddr',
                             'omega3_noala_ffq', 'epa_ffq', 'dpa_ffq', 'dha_ffq',
                             'omega3_noala10v', 'epa10v', 'dpa10v', 'dha10v'
                             ) ])
colSums(missing1)
sum(rowSums(missing1)!=0)



##### metadata by stool ####


write.table(metadata_intrst, file = "./data/meta_stn.pcl",  sep = '\t', quote = F, eol = '\n',  row.names=F)
write.table(logmeta, file = "./data/logmeta_stn.pcl",  sep = '\t', quote = F, eol = '\n',  row.names=F)
write.table(totmeta, file = "./data/totmeta_stn.pcl",  sep = '\t', quote = F, eol = '\n',  row.names=F)

write.table(logspecies, file = "./data/logspecies_filt.pcl",  sep = '\t', quote = F, eol = '\n',  row.names=T)
write.table(bispecies, file = "./data/bispecies_filt.pcl",  sep = '\t', quote = F, eol = '\n',  row.names=T)

########## DICTIONARY #############

##How to import?
species <-read.table(file = './data/species_filt.csv',
                     sep = ',',  row.names = 1,  header = TRUE,   check.names = FALSE)
logspecies <-read.table(file = './data/logspecies_filt.pcl', row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
bispecies <- read.table('./data/bispecies_filt.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")

meta_stn <- read.table('./data/meta_stn.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")
logmeta_stn <- read.table('./data/logmeta_stn.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")



reg <- lmer(data = metadata_intrst, logcrp_plasma ~ omega3_noala_ffq + (1|ID1))


anova(reg)

