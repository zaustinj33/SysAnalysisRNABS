library(dplyr)
library(tidyr)
library(ggpubr)
##%

## Load counts
raw_counts <- read.csv("featureCounts_rawCountsMT_clean.txt",
                       header = T, sep = '\t', stringsAsFactors = F, row.names = 1)
raw_counts[is.na(raw_counts)] <- 1
raw_counts[raw_counts == 0] <- 1

# Create normalized counts table
dds.cpm <- apply(raw_counts,2, function(x) ((x/sum(x))*1000000))
dds.cpm <- dds.cpm[,c(5,1,2,3,4)]
cpm.long <- pivot_longer(as.data.frame(dds.cpm), cols = everything())
# Boxplot of CPM (counts per million)
ggboxplot(cpm.long, x='name', y='value', xlab="", ylab="Counts per million",las=2,
         fill = 'name', palette = c('blue', 'orange', 'green', 'red', 'purple')) +
  yscale('log10', .format = TRUE)

abline(h=median(dds.cpm),col="blue")
#title("Boxplots of logCPMs")
ggsave(filename = "CPM_barplot.png",height = 4,width = 6,dpi = 200)

# Scatter plot comparing CPM of conditions
conditions <- c('G1','G2','G3','G4')
compairson_plot <- function(condition){
  plot_df <- as.data.frame(dds.cpm) %>% select(condition, 'G5')
  plot_df$lFC <- ifelse((log2(plot_df[[1]]) - log2(plot_df[[2]])) > 2,
  'Over', ifelse((log2(plot_df[[1]]) - log2(plot_df[[2]])) < -2, 'Under', 'None'))
  print(table(plot_df$lFC))

  ggscatter(as.data.frame(plot_df), y=colnames(plot_df)[[1]], x='G5', color = 'lFC', palette = c("grey", 'blue', 'orange'),
            conf.int = TRUE, cor.coef = TRUE) +
    xscale('log10', .format = TRUE) +
    yscale('log10', .format = TRUE)

  ggsave(filename = paste0("lFC_cpm_",condition, ".png") ,height = 4,width = 4,dpi = 400)

}
lapply(conditions, compairson_plot)
