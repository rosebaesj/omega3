options(max.print=1000000)
library(data.table)
library(readr)
library(haven)
library(plyr)
library(dplyr)
library(vegan)
library(scales)
library(RColorBrewer)
library(grid)
library(pheatmap)
library(lme4)
library(ggplot2)
library(cowplot)
library(ggsignif)
library(car)
library(Maaslin2) #Maaslin2     * 0.99.12     2020-06-18
library(gridExtra)
library(knitr)
library(readxl)
library(ggpubr)
#update.packages()


#functions i will need later 
## Edit body of pheatmap:::draw_colnames, customizing it to your liking
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}

## For pheatmap_1.0.8 and later:
## I change rot=45 to rot=90
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))


maaslin_heatmap <- function(
  maaslin_output, 
  title = "", cell_value = "qval", data_label = 'data', metadata_label = 'metadata', 
  gaps_row = gaps, gaps_col = gapscol, border_color = "black", 
  paletteLength= 50,
  colors = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(paletteLength), annotation_row, ...) #)
{
  # read MaAsLin output
  df <- read.table( 
    maaslin_output,
    header = TRUE, sep = "\t", fill = TRUE, comment.char = "" , check.names = FALSE) 
  metadata <- df$metadata
  data <- df$feature
  value <- NA
  # values to use for coloring the heatmap
  if (cell_value == "pval"){
    value <- -log(df$pval)*sign(df$coef)
    value <- pmax(-4, pmin(4, value))
  }else if(cell_value == "qval"){
    value <- -log(df$qval)*sign(df$coef)
    value <- pmax(-6, pmin(6, value))
  }else if(cell_value == "coef"){
    value <- df$coef
    value <- pmax(-0.01, pmin(0.01, value))
  }
  n <- length(unique(metadata))
  m <- length(unique(data))
  a = matrix(0, nrow=m, ncol=n)
  for (i in 1:length(data)){
    # if (abs(a[as.numeric(metadata)[i], as.numeric(metadata)[i]]) > abs(value[i]))
    # next
    a[as.numeric(data)[i], as.numeric(metadata)[i]] <- value[i]
  }
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks <- c(seq(-abs(max(min(value),max(value))), 0, length.out=ceiling(paletteLength/2) +1),
                seq(abs(max(min(value),max(value)))/paletteLength, abs(max(min(value),max(value))), length.out=floor(paletteLength/2)))
  # noasterisk <- gsub('*', '', levels(data))  
  # rownames(a) <- substr(levels(data), 0, nchar(noasterisk))
  
  rownames(a) <- substr(levels(data), 0, nchar(levels(data)))
  colnames(a) <- substr(levels(metadata), 3, nchar(levels(metadata))) # to maintain order by numbering# 
  #colnames(a) <- levels(metadata) # if just row names #
  
  #rownames(annotation_row) <- rownames(a)
  
  #colnames(a) <-  sapply(colnames(a), pcl.sub )
  
  p <- pheatmap(a, cellwidth = 14, cellheight = 14,   # changed to 3
                main = title,
                fontsize = 10,
                kmeans_k = NA,
                border=TRUE,
                show_rownames = need_rownames, show_colnames = T,
                scale="none",
                clustering_method = "complete",
                cluster_rows = FALSE, cluster_cols = FALSE,
                # clustering_distance_rows = "euclidean", 
                # clustering_distance_cols = "euclidean",
                legend=needlegend,
                border_color = border_color,
                color = colors,
                breaks=myBreaks,
                treeheight_row=0,
                treeheight_col=0,
                gaps_row=gaps,
                gaps_col=gapscol,
                # annotation_row=annotation_row,
                annotation_legend=annoleg,
                # cutree_rows = clustnum,
                na_col='#F4F6F6',
                fontsize_number=15,
                # display_numbers = matrix(ifelse(a < 0.5, '*', ifelse(a < 0.0, '', '')),  nrow(a)), ...)
                # display_numbers = matrix(ifelse(a > 0, "+", ifelse(a < 0.0, "-", "")),  nrow(a)))
                display_numbers = matrix(ifelse(df$qval < 0.25, "*", ifelse(df$qval > 0.25, "", "")), nrow(a)))
  
  return(p)
}

nature_theme <- theme(axis.text.x = element_text(size = 8, vjust = 1),
                      axis.text.y = element_text(size = 8, hjust = 1),
                      axis.title=element_text(size = 10  ),
                      plot.title =element_text(size=7, face='bold'),
                      legend.title=element_text(size=6, face='bold'),
                      legend.text=element_text(size=6),
                      axis.line = element_line(colour = 'black', size = .25),
                      axis.line.x = element_line(colour = 'black', size = .25), axis.line.y = element_line(colour = 'black', size = .25)) 

# reassign the key for all 4 species files 
reassignKey = function(myDataName)
{
  newDataName <- merge(myDataName, z_idkey,  by="id")
  newDataName <- newDataName[order(newDataName$ID1), ]
  rownames(newDataName) <- newDataName$ID1
  newDataName <- newDataName[ , -which(names(newDataName) %in% c('id', 'ID1'))] 
  return(newDataName)    
}
