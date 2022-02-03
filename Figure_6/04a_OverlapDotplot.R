library(ggpubr)
library(ggplot2)
library(MASS)
library(ggsci)
library(dplyr)

## Dotplot showing similarity of methylation level between two replicates for 
# all transcripts
# set directories if necessary #
# parent directory
setwd(paste0(getwd(), "/Chapter_5"))
save_dir <- paste0(getwd(),"/figures")

#Overlap from comprehensive union file
all_m5C <- read.csv("allm5C_libraries_filteredDepth.csv")

## MANUALLY CHANGE SAVE FILE NAME AND MT FILTER IN FUNCTION ##
conditions <- c("G","MF","pMF","SRR","")

# color density function
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

## dotplot ##
overlap_files <- list.files(pattern = "_allOverlapDepth.csv", full.names = T)
#View(read.csv(overlap_files[[1]]))

## Comparing replicates
replicate_dotplot <- function (overlap) {
  name <- gsub("_allOverlapDepth.csv","",basename(overlap))
  all_overlap <- read.csv(overlap)

  overlap_df <- all_overlap[,grepl(paste0(name), colnames(all_overlap))]
  # only grab methylLevel columns
  methylDF <- data.frame(group=all_overlap$group,
                         ML_1=as.numeric(overlap_df[[3]]),
                         ML_2=as.numeric(overlap_df[[6]])
  )
  methylDF <- methylDF[(methylDF[,2] >= 0.1) | (methylDF[,3] >= 0.1),]
  methylDF <- methylDF[complete.cases(methylDF),]
  print(nrow(methylDF))
  methylDF$density <- get_density(methylDF$ML_1, methylDF$ML_2, n = nrow(methylDF))

  p<-ggplot(methylDF) + 
    geom_point(aes(x = ML_1, y = ML_2, color = density),size = 3) + theme_bw() +
    geom_abline(slope = 1, linetype = 'longdash', size = 1.5) +
    scale_color_gradientn(colors = c('blue','yellow','#DC0000FF'), values = c(0,0.5,1)) +
    theme(legend.title = element_blank(), legend.position = 'none', axis.text = element_text(size = 18),
          text = element_text(size = 18, colour = 'black'), axis.line = element_line(size = 1),
          panel.grid=element_blank(), panel.border = element_rect(size = 2), axis.ticks.length=unit(.25, "cm"),
          plot.margin = margin(1,1,1,1,'cm')) +
    scale_x_continuous(expand = c(0,0), breaks = c(0,0.5,1),limits = c(0:1)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0,0.5,1),limits = c(0:1)) +
    #ggtitle(name) +
    xlab("m5c level of rep1") + 
    ylab("m5c level of rep2") +
    stat_cor(method="pearson",aes(x = ML_1, y = ML_2),label.x.npc = 0.05)
  p
  ggsave(path = save_dir,filename = paste0("/Dotplot_", name, "_noFDRdepth.png")
    ,plot = p,height = 5,width = 5,dpi = 400)
}
plot <- replicate_dotplot(overlap_files[1])

lapply(overlap_files, replicate_dotplot)

#%%

## dotplot ## 20x only filter MT libs
overlap_files <- list.files(pattern = "G_allUnionDepth.csv", full.names = T)
#View(read.csv(overlap_files[[1]]))

conditions <- c("G1|G2", "G1|G3", "G1|G4", "G2|G3", "G2|G4", "G3|G4")
preDMS_all <- read.csv("allm5C_libraries.csv", comment.char = '',header = T)
preDMS_all[,c(4:10,14:19,23:28,32:37)] <- NULL
preDMS_all[,4] <- NULL
## Comparing replicates
replicate_dotplot <- function (condition) {
  name <- gsub("\\|","_", condition)
  print(name)
  all_overlap <- read.csv(overlap_files)

  # Comparison condition
  overlap_df <- preDMS_all[,grepl(condition, colnames(all_overlap))]
  overlap_df <- overlap_df[(overlap_df[,1] >=20) & (overlap_df[,4] >=20) &
                             (overlap_df[,2] >=1) & (overlap_df[,5] >=1),]
  overlap_df <- overlap_df[complete.cases(overlap_df),]
  write.csv(overlap_df, paste0("Dotplot_", name, "_compare.csv"))
  # only grab methylLevel columns
  methylDF <- data.frame(
                         ML_1=as.numeric(overlap_df[[3]]),
                         ML_2=as.numeric(overlap_df[[6]])
  )
  methylDF <- methylDF[complete.cases(methylDF),]
  print(nrow(methylDF))
  #methylDF$density <- get_density(methylDF$ML_1, methylDF$ML_2, n = nrow(methylDF))

  p<-ggplot(methylDF) +
    geom_point(aes(x = ML_1, y = ML_2),size = 1) + theme_bw() +
    geom_abline(slope = 1, linetype = 'longdash', size = 1.5) +
    scale_color_gradientn(colors = c('blue','yellow','#DC0000FF'), values = c(0,0.5,1)) +
    theme(legend.title = element_blank(), legend.position = 'none', axis.text = element_text(size = 18),
          text = element_text(size = 18, colour = 'black'), axis.line = element_line(size = 1),
          panel.grid=element_blank(), panel.border = element_rect(size = 2), axis.ticks.length=unit(.25, "cm"),
          plot.margin = margin(1,1,1,1,'cm')) +
    scale_x_continuous(expand = c(0,0), breaks = c(0,0.5,1),limits = c(0:1)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0,0.5,1),limits = c(0:1)) +
    ggtitle(condition) +
    xlab("m5c level of rep1") +
    ylab("m5c level of rep2")
    #stat_cor(method="pearson",aes(x = ML_1, y = ML_2),label.x.npc = 0.05)
  p
  #ggsave(path = save_dir,filename = paste0("/Dotplot_", name, "_compare.png")
   # ,plot = p,height = 5,width = 5,dpi = 500)
}
lapply(conditions, replicate_dotplot)

sample <- "G1"
preDMS_MT <- preDMS_all[preDMS_all$chrom == 'MT',]
overlap_df <- preDMS_MT[,grepl(sample, colnames(all_overlap))]
overlap_df <- overlap_df[complete.cases(overlap_df),]
overlap_df$sample <- sample
overlap_df$bin_cov <- ifelse(overlap_df$cov_G1 > 10000, "1) > 10,000",
                             ifelse(overlap_df$cov_G1 > 1000, "2) 10,000 > x > 1,000",
                             ifelse(overlap_df$cov_G1 > 100, "3) 1,000 > x > 100",
                             ifelse(overlap_df$cov_G1 > 20,  "4) < 100","5) < 20"))))

ggplot(aes(y=methRate_G1, x=sample, fill=bin_cov), data=overlap_df) + geom_boxplot()
ggsave(path = save_dir,filename = "G1_MT_MLbins.png" ,height = 5,width = 5,dpi = 500)
