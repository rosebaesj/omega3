# [2] metadata set
### (1) cohort FFQ cumulative average
'ala10v','epa10v','dha10v', 'dpa10v', 'trans10v', 'omega610v', 'omega310v', 'omega3_noala10v'
Note: I tried further including act10 and bristol score as continuous in the model, which did not influence the results too much.
~~~~~~~~~~~~r~~~~~~~~~~~~  
ffqcum_keep <- c('participant', 'agemlvs', 'calor10v', 'abx_avg')
ffqcum_keep_nocal <- c('participant', 'agemlvs', 'abx_avg')

ala_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep, 'ala10v')] 
ala_nocal_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep_nocal, 'ala10v')] 
epa_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep, 'epa10v')] 
epa_nocal_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep_nocal, 'epa10v')] 
dha_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep, 'dha10v')] 
dha_nocal_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep_nocal, 'dha10v')] 
dpa_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep, 'dpa10v')] 
dpa_nocal_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep_nocal, 'dpa10v')] 
trans_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep, 'trans0v')] 
trans_nocal_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep_nocal, 'trans10v')] 
omega6_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep, 'omega610v')] 
omega6_nocal_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep_nocal, 'omega610v')] 
omega3_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep, 'omega310v')] 
omega3_nocal_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep_nocal, 'omega310v')] 
omega3_noala_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep, 'omega3_noala10v')] 
omega3_noala_nocal_ffqcum <- z_merged.all[ , colnames(z_merged.all) %in% c(ffqcum_keep_nocal, 'omega3_noala10v')] 
~~~~~~~~~~~~~~~~~~~~~~~~  
### (2) MLVS FFQ
'ala_avg', 'epa_avg', 'dha_avg', 'dpa_avg', 'omega3_avg', 'omega3_noala_avg', 'omega6_avg','trans_avg'
~~~~~~~~~~~r~~~~~~~~~~~~~  
ffq_keep <- c('participant', 'agemlvs', 'calor_avg', 'abx_avg')
ffq_keep_nocal <- c('participant', 'agemlvs', 'abx_avg')

ala_avg_fs_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep, 'ala_avg')] 
ala_avg_fs_nocal_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep_nocal, 'ala_avg')] 
epa_avg_fs_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep, 'epa_avg')] 
epa_avg_fs_nocal_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep_nocal, 'epa_avg')] 
dha_avg_fs_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep, 'dha_avg')] 
dha_avg_fs_nocal_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep_nocal, 'dha_avg')] 
dpa_avg_fs_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep, 'dpa_avg')] 
dpa_avg_fs_nocal_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep_nocal, 'dpa_avg')] 
omega3_avg_fs_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep, 'omega3_avg')] 
omega3_avg_fs_nocal_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep_nocal, 'omega3_avg')] 
omega3_noala_avg_fs_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep, 'omega3_noala_avg')] 
omega3_noala_avg_fs_nocal_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep_nocal, 'omega3_noala_avg')] 
omega6_avg_fs_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep, 'omega6_avg')] 
omega6_avg_fs_nocal_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep_nocal, 'omega6_avg')] 
trans_avg_fs_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep, 'trans_avg')] 
trans_avg_fs_nocal_ffq <- z_merged.all[ , colnames(z_merged.all) %in% c(ffq_keep_nocal, 'trans_avg')] 
~~~~~~~~~~~~~~~~~~~~~~~~ 

### (3) MLVS DDR
#'a_omega3_fs', 'a_omega3_noala_fs', 'a_ala_fs_dr_w1', 'a_f205_fs_dr_w1', 'a_f226_fs_dr_w1', 'a_p22_5_fs_dr_w1', 'a_trn07_fo_dr_w1'
~~~~~~~~~~~r~~~~~~~~~~~~~ 
ddr_keep <- c('participant', 'agemlvs', 'calor_fs_dr', 'abx')
ddr_keep_nocal <- c('participant', 'agemlvs', 'abx')
a_omega3_fs_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_omega3_fs_dr')] 
a_omega3_fs_nocal_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep_nocal, 'a_omega3_fs_dr')] 
a_omega3_noala_fs_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_omega3_noala_fs_dr')] 
a_omega3_noala_fs_nocal_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep_nocal, 'a_omega3_noala_fs_dr')] 
a_ala_fs_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_ala_fs_dr')] 
a_ala_fs_nocal_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep_nocal, 'a_ala_fs_dr')] 
a_trans_fs_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_trn07_fo_dr')] 
a_trans_fs_nocal_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep_nocal, 'a_trn07_fo_dr')] 
a_epa_fs_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_f205_fs_dr')] 
a_epa_fs_nocal_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep_nocal, 'a_f205_fs_dr')] 
a_dha_fs_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_f226_fs_dr')] 
a_dha_fs_nocal_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep_nocal, 'a_f226_fs_dr')] 
a_dpa_fs_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_p22_5_fs_dr')] 
a_dpa_fs_nocal_ddr  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep_nocal, 'a_p22_5_fs_dr')] 
omega3_tfat_ddr <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_omega3_fs_dr', "a_fat_fs_dr")] 
trans_tfat_ddr <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_trans_fs_dr', "a_fat_fs_dr")] 
epa_tfat_ddr <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_epa_fs_dr', "a_fat_fs_dr")] 
dha_tfat_ddr <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_dha_fs_dr', "a_fat_fs_dr")] 
dpa_tfat_ddr <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(ddr_keep, 'a_dpa_fs_dr', "a_fat_fs_dr")] 
~~~~~~~~~~~~~~~~~~~~~~~~ 

### (4) MLVS biomarker
~~~~~~~~~~~r~~~~~~~~~~~~~ 
biomarker_keep <- c('participant', 'agemlvs', 'abx', 'bmi12')
biomarker_logcrp  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logcrp')] 
biomarker_logcrp_complete <- biomarker_logcrp[!is.na(biomarker_logcrp$logcrp), ]
biomarker_loghdl  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'loghdl')] 
biomarker_loghdl_complete <- biomarker_loghdl[!is.na(biomarker_loghdl$loghdl), ]
biomarker_logtg  <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtg')] 
biomarker_logtg_complete <- biomarker_logtg[!is.na(biomarker_logtg$logtg), ]

a_omega3_fs_ddr_logcrp <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logcrp','a_omega3_fs_dr')] 
a_omega3_fs_ddr_logtchdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtchdl','a_omega3_fs_dr')] 
a_omega3_fs_ddr_loghdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'loghdl','a_omega3_fs_dr')] 
a_omega3_fs_ddr_logtg <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtg','a_omega3_fs_dr')] 

a_omega3_noala_fs_ddr_logcrp <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logcrp','a_omega3_noala_fs_dr')] 
a_omega3_noala_fs_ddr_logtchdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtchdl','a_omega3_noala_fs_dr')] 
a_omega3_noala_fs_ddr_loghdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'loghdl','a_omega3_noala_fs_dr')] 
a_omega3_noala_fs_ddr_logtg <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtg','a_omega3_noala_fs_dr')] 

a_ala_fs_ddr_logcrp <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logcrp','a_ala_fs_dr')] 
a_ala_fs_ddr_logtchdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtchdl','a_ala_fs_dr')] 
a_ala_fs_ddr_loghdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'loghdl','a_ala_fs_dr')] 
a_ala_fs_ddr_logtg <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtg','a_ala_fs_dr')] 

a_trans_fs_ddr_logcrp <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logcrp','a_trans_fs_dr')] 
a_trans_fs_ddr_logtchdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtchdl','a_trans_fs_dr')] 
a_trans_fs_ddr_loghdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'loghdl','a_trans_fs_dr')] 
a_trans_fs_ddr_logtg <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtg','a_trans_fs_dr')] 

a_epa_fs_ddr_logcrp <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logcrp','a_epa_fs_dr')] 
a_epa_fs_ddr_logtchdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtchdl','a_epa_fs_dr')] 
a_epa_fs_ddr_loghdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'loghdl','a_epa_fs_dr')] 
a_epa_fs_ddr_logtg <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtg','a_epa_fs_dr')] 

a_dha_fs_ddr_logcrp <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logcrp','a_dha_fs_dr')] 
a_dha_fs_ddr_logtchdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtchdl','a_dha_fs_dr')] 
a_dha_fs_ddr_loghdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'loghdl','a_dha_fs_dr')] 
a_dha_fs_ddr_logtg <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtg','a_dha_fs_dr')] 

a_dpa_fs_ddr_logcrp <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logcrp','a_dpa_fs_dr')] 
a_dpa_fs_ddr_logtchdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtchdl','a_dpa_fs_dr')] 
a_dpa_fs_ddr_loghdl <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'loghdl','a_dpa_fs_dr')] 
a_dpa_fs_ddr_logtg <- z_w12.all_species.all[ , colnames(z_w12.all_species.all) %in% c(biomarker_keep, 'logtg','a_dpa_fs_dr')] 
~~~~~~~~~~~~~~~~~~~~~~~~ 


# drop certain values
## 1) drop observations with missing short-term fiber
I don't know why I need to use the complete dataset to run MaAsLin2 in my mac...
```diff
- what dose this mean?
```
~~~~~~~~~~~r~~~~~~~~~~~~~ 
a_omega3_fs_ddr_complete <- a_omega3_fs_ddr[a_omega3_fs_ddr$a_omega3_fs_dr>0,]
a_omega3_fs_ddr_complete <- na.omit(a_omega3_fs_ddr_complete)
~~~~~~~~~~~~~~~~~~~~~~~~ 

## 2) Tobyn
```diff
- what dose this mean?
```
#### check distribution of fiber 
~~~~~~~~~~~~r~~~~~~~~~~~~
graphics.off()
histogram(a_omega3_fs_ddr$a_omega3_fs_dr)
quantile(a_omega3_fs_ddr$a_omega3_fs_dr,  probs = seq (0,1,0.1), na.rm = TRUE,
         names = TRUE)
~~~~~~~~~~~~~~~~~~~~~~~~
0%       10%       20%       30%       40%       50%       60%       70%       80%       90%      100% 
9.392454 16.591799 19.089333 20.941394 22.565379 24.031779 25.967835 27.866991 30.168380 35.819875 55.518952 


#### IQR: Compute the interquartile range for a vector.
~~~~~~~~~~~~r~~~~~~~~~~~~
IQR(a_omega3_fs_ddr$a_omega3_fs_dr, na.rm = TRUE) 

quantile(a_omega3_fs_ddr$a_omega3_fs_dr,  0.75, na.rm = TRUE,
         names = TRUE)
~~~~~~~~~~~~~~~~~~~~~~~~
1.085002

75% 
2.75162

outlier: Q3+1.5IQR=2.75162+1.5*1.085002=4.379123
Q3+3IQR=2.75162+3*1.085002=6.006626

#### drop observations with short-term fiber >Q3+1.5IQR
a_omega3_fs_ddr_1.5IQR <- a_omega3_fs_ddr[a_omega3_fs_ddr$a_omega3_fs_dr<=4.379123,]
a_omega3_fs_ddr_3IQR <- a_omega3_fs_ddr[a_omega3_fs_ddr$a_omega3_fs_fs_dr<=6.006626,]

## what is this
~~~~~~~~~~~~r~~~~~~~~~~~~
mlvs_metadata_925 <- z_w12.all_species.all[,!colnames(z_w12.all_species.all) %in% colnames(species_all)]
~~~~~~~~~~~~~~~~~~~~~~~~
```diff
- what dose this mean?
```
