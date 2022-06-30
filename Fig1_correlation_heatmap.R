

#### Correlation analysis ####

cor_data <- subset(meta_species, 
                   select=c(agemlvs, bmi12, calor10v, logcrp_w1, loghdl_w1, logtg_w1))
cor <- rcorr(as.matrix(cor_data ),type=c("spearman"))
cor

meta_species$'Age' <- meta_species$agemlvs
meta_species$'BMI' <- meta_species$bmi12
meta_species$'Calory_intake' <- meta_species$calor10v
meta_species$'Plasma_CRP' <- meta_species$logcrp_w1
meta_species$'Plasma_HDL' <- meta_species$loghdl_w1
meta_species$'Plasma_TG' <- meta_species$logtg_w1

list <- c('Age', 
          'BMI', 
          'Calory_intake', 
          'Plasma_CRP', 
          'Plasma_HDL', 
          'Plasma_TG')

pawt_data <- subset(meta_species, select=list)

cormat <- rcorr(as.matrix(pawt_data),type=c("spearman"))
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
# Heatmap

level_x_order <- factor(melted_cormat$variable, level = rev(list))

level_y_order <- factor(melted_cormat$var1, level = list)

ggheatmap <- ggplot(melted_cormat, aes(level_x_order, level_y_order, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  geom_text(aes(label=value), color="black", size=4, show.legend = TRUE) +
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




