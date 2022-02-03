## Generates sites list for all transcripts in preDMS table, or all potential important
## sites. Take all_preDMS table as input, and outputs an annotated table with site
## description (3UTR, 5UTR, or exon), and bin location as input for plot_site_bins.R 

library(ensembldb)
library(EnsDb.Mmusculus.v79)
library(data.table)
library(dplyr)

setwd(paste0(getwd(),"/Chapter_5"))
home_dir <- getwd()
save_dir <- getwd()

## Create site list annotations for all sites contained in preDMS table
# WARNING: this is a computationally heavy step and will eat upwards of 16GB of RAM
# if you must use it, plan to break input into small chunks (5000 rows each) and merge
# at end of step

# import filtered m5c list
preDMS_all <- read.csv("allm5C_libraries_filteredDepth.csv", comment.char = '',header = T)
#preDMS_all$geno_group  <- paste(preDMS_all$, preDMS_all$transcript_start)
# Create test dataframe
#preDMS_test <- preDMS_all[sample(nrow(preDMS_all), 100), ]

# import reference
edbx <- EnsDb.Mmusculus.v79
prep_annotations <- function (preDMS_file){
  # read in Input files
  name <- "all_merged"
  overlap <- preDMS_file
  overlap_genome <- data.frame(chrom=overlap$chrom, position=overlap$refPos, group=overlap$group)
  #annotation

  dat1 = overlap_genome
  dat2 <- read.delim("mm10_anno_79_genes.txt")
  
  final_peak_diff_list = seq(1, dim(dat1)[1])
  common_peaks = list()
  tmp = c(); j=1
  
  ## calculate set diff
  for(ipeak in 1:dim(dat2)[1]){
    if(dim(dat1)[1]>1){
      indx = which(dat1$chrom %in% dat2$chrom[ipeak] & (dat1$position >= dat2$start[ipeak])& (dat1$position <= dat2$end[ipeak])) 
      if(length(indx)>0){
        common_peaks [[j]] = cbind.data.frame(dat1[indx,],dat2[ipeak,])
        j = j+1
        dat1 = dat1[ -indx,]
        rownames(dat1) = NULL
      }
    }
  }
  tmp[4] <- NULL
  tmp = do.call('rbind', common_peaks)
  
  
  #create columns for bin-mapping 
  #length of transcript region: (stop-start), 
  tmp$length <- tmp$end - tmp$start
  # bin size: 5'UTR=5,CDS=22,3'UTR=18; numbers derived from Jialin code
  # create filler term
  tmp$binSize <- 0
  tmp$binSize <- ifelse(tmp$type == "3UTR", ((tmp$length-1)/18), tmp$binSize)
  tmp$binSize <- ifelse(tmp$type == "5UTR", ((tmp$length-1)/5), tmp$binSize)
  tmp$binSize <- ifelse(tmp$type == "exon", ((tmp$length-1)/22), tmp$binSize)
  
  #bin location: (chrom. position - (start +1) / [bin size])
  tmp$binLoc <- ceiling((tmp$position - (tmp$start)+1)/tmp$binSize)
  
  anno <- tmp
  #write.csv(anno, "all_anno.csv",row.names = F)
  write.csv(dat1, "all_leftover.csv",row.names = F)
  
  return(anno)
}

tmp_fxnTest <- prep_annotations(preDMS_all)
tmp_fxnTest[5] <- NULL # remove extra chrom column
#write.csv(tmp_fxnTest, paste0(save_dir,name,"_m5c_genome_ANNO.csv"))

all_merged_anno <- right_join(preDMS_all, tmp_fxnTest, by='group')
all_merged_anno$refStrand <- NULL
write.csv(all_merged_anno,'allm5C_libraries_filteredDepthAnno.csv')

