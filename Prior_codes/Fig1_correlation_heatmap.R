species <-read.table(file = './data/species_filt.csv',
                     sep = ',',  row.names = 1,  header = TRUE,   check.names = FALSE)

meta_stn <- read.table('./data/meta_stn.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")

totmeta_stn <- read.table('./data/totmeta_stn.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")


logspecies <-read.table(file = './data/logspecies_filt.pcl', row.names = 1,  header = TRUE, sep='\t',  check.names = FALSE, quote ="")
bispecies <- read.table('./data/bispecies_filt.pcl',row.names = 1, header=TRUE, sep='\t', check.names=TRUE, quote ="")


meta_logspecies <- merge(totmeta_stn, logspecies, by = 0)
meta_bispecies <- merge(totmeta_stn, bispecies, by = 0)

#### Correlation analysis ####

cor_data <- subset(meta_species, 
                   select=c(agemlvs, bmi12, calor10v, logcrp_w1, loghdl_w1, logtg_w1))
cor <- rcorr(as.matrix(cor_data ),type=c("spearman"))
cor

list <- c('agemlvs', 'bmi12', 'logcrp_plasma','logAA_pfa',
          'omega3_ddr', 
          'omega3_ffq',
          'omega310v',
          #'omega3_pfa',
          'omega3_noala_ddr',# 'epa_ddr', 'dpa_ddr', 'dha_ddr', 'ala_ddr', 'omega3_ddr',
          'omega3_noala_ffq',#'epa_ffq', 'dpa_ffq', 'dha_ffq','ala_ffq','omega3_ffq',
          'omega3_noala10v',#'epa10v', 'dpa10v', 'dha10v','ala10v','omega310v',
          'omega3_noala_pfa','epa_pfa', 'dpa_pfa', 'dha_pfa', 'ala_pfa','omega3_pfa',
          'omega6_pfa','sat_pfa', 'monounsat_pfa', 'trans_pfa'
)
list<- c('bmi_paq1', 'wt_paq1', 'waist_paq1', 'whr_paq1', 'wtchg',
         'logcrp_plasma', 'hba1c')
corrmat <- totmeta_stn[,c('agemlvs', 'bmi12', 
                          'logcrp_plasma','logAA_pfa',
                          'omega3_ddr',
                          'omega3_ffq',
                          'omega310v',
                          'omega3_pfa',
                          'omega3_noala_ddr', #'epa_ddr', 'dpa_ddr', 'dha_ddr', 'ala_ddr', 'omega3_ddr',
                          'omega3_noala_ffq',#'epa_ffq', 'dpa_ffq', 'dha_ffq','ala_ffq','omega3_ffq',
                          'omega3_noala10v',#'epa10v', 'dpa10v', 'dha10v','ala10v','omega310v',
                          'omega3_noala_pfa')]#,'epa_pfa', 'dpa_pfa', 'dha_pfa', 'ala_pfa','omega3_pfa')]


corrmat <- totmeta_stn[,c('agemlvs', 'bmi12', 
                          'logcrp_plasma','logAA_pfa',
                          'omega3_pfa','omega6_pfa', 
                          'sat_pfa', 'monounsat_pfa', 'trans_pfa')]




cormat <- rcorr(as.matrix(corrmat),type=c("spearman"))

cormat
rcx <- cormat
str(rcx)
df.rcx<-data.frame(rcx$r)
df.rcx<-round(df.rcx,2)
df.rcx

# Get lower triangle of the correlation matrix
get_lower_tri<-function(df.rcx){
  df.rcx[lower.tri(df.rcx)] <- NA
  return(df.rcx)
}
lower_tri <- get_lower_tri(df.rcx)
lower_tri
lower_tri$var1<-rownames(lower_tri)

melted_cormat <- melt(lower_tri, id.vars='var1', na.rm = TRUE)
melted_cormat




pcx <- cormat
str(pcx)
df.pcx<-data.frame(pcx$P)
df.pcx<-round(df.pcx,4)
df.pcx

# Get lower triangle of the correlation matrix
get_lower_tri<-function(df.rcx){
  df.rcx[lower.tri(df.rcx)] <- NA
  return(df.rcx)
}
lower_tri_p <- get_lower_tri(df.pcx)
lower_tri_p
lower_tri_p$var1<-rownames(lower_tri_p)


melted_cormat_p <- melt(lower_tri_p, id.vars='var1', na.rm = TRUE)
melted_cormat_p$stars<-cut(as.numeric(melted_cormat_p$value), breaks=c(-Inf, 0.005, 0.01, 0.05, Inf), label=c("***","**", "*", ""))


melted_rp<- left_join(melted_cormat, melted_cormat_p, by = c("var1", "variable"))

# Heatmap

level_x_order <- factor(melted_cormat$variable, level = rev(list))

level_y_order <- factor(melted_cormat$var1, level = list)

ggplot(melted_rp, aes(level_x_order, level_y_order, fill = value.x))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  geom_text(aes(label=stars), color="black", size=4, show.legend = TRUE) +
  theme_minimal()+ # minimal theme
  theme(legend.title = element_text(size = 10,color="black"),
        legend.text = element_text(size = 10,color="black"),
        axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1,color="black"),
        axis.text.y = element_text(size = 10,color="black"),
        axis.title.x= element_blank(),
        axis.title.y= element_blank())+
  coord_fixed()
# Print the heatmap
print(ggheatmap)


png(file="./figure_generated/fig1_correlation_pawt.png",
    width=5,height=5, units="in", res = 1000)
grid.arrange(ggheatmap, ncol=1, nrow=1)
dev.off()




