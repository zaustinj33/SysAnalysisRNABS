library(ggseqlogo)
library(ggplot2)
library(ggsci)
library(dplyr)
###################################

# Import overlap files
setwd(paste0(getwd(), "/Chapter_5"))
save_dir <- paste0(getwd(),"/figures")

#Overlap from coprehensive union file
all_m5C <- read.csv("allm5C_libraries_filteredDepthAnno.csv")
## MANUALLY CHANGE SAVE FILE NAME AND MT FILTER IN FUNCTION ##
conditions <- c("G","MF","pMF","SRR","")

#overlap_files <- list.files(pattern = "overlap_filtered.csv$", recursive = T) # Raw counts
overlap_files <- list.files(pattern = "_allOverlapDepth.csv", full.names = T)  # Overlap
#overlap_files <- list.files(pattern = "allUnion.csv", full.names = T)  # Union


# Define colour scheme
cs1 <- make_col_scheme(chars=c('A', 'U', 'C', 'G'),
                      #cols=c(pal_npg("nrc")(10)[c(4,2,8)],"#fbc02d"))
                      cols=c(pal_npg("nrc")(10)[c(1,8)],"#0000FF","#fbc02d"))

#WT <- read.csv(overlap_files[[3]])

create_seqplot <- function (overlap_file) {
  overlap <- read.csv(union_files[[2]])
  name <-  gsub("overlap.*$","",basename(overlap_file))
  print(paste('number of sites in ', name, ': ',nrow(overlap)))
  
  #MT filter if needed
  #overlap <- overlap[overlap$group %in% MT_list,]
  
  # Create Seq df
  seqLogoDF <- data.frame(group = overlap$group, seqContext = as.character(overlap$seqContext))
  # get lengths of seqlogo for some reason
  seqLogoDF$length <- nchar(as.character(seqLogoDF[,1])) ##defining sequence context length
  #a <- seqLogoDF[seqLogoDF$length == 21,] ##only choose context with 21 bp
  #substring(seqLogoDF$seqContext.x,11) <- 'mC'
  #b <- toupper(a[,1]) ##make lowercase uppercase
  b <- toupper(seqLogoDF$seqContext)
  b <- gsub("T", "U", b)
  
  g <- ggseqlogo(b, seq_type="rna", method = "prob",col_scheme = cs1) +
    theme_bw() + xlab('Position') + ggtitle('Huang') +
    theme(legend.title = element_blank(), legend.position = 'none', 
          axis.text = element_text(size = 8), text = element_text(size = 18, colour = 'black'),
          axis.line = element_line(size = 1),
          panel.grid=element_blank(), panel.border = element_blank(), axis.ticks.length=unit(.25, "cm"),
          plot.margin = margin(1,1,1,1,'cm'))
  g
  
  
  ggsave(path = save_dir,filename =  paste0("SeqContext_",
                                            'Huang', "_HighConf.png"),plot = g,height = 3,width = 5,dpi = 500)
}
create_seqplot(overlap_files[[1]])
lapply(overlap_files[c(2,3,4)], create_seqplot)
