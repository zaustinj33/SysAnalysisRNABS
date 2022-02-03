library('devtools')
library('GenomicRanges')
library('biomaRt')
library('GenomicFeatures')
library(dplyr)
#%%
## rBS GRanges format
#seqnames: chrX
#ranges: chromosome position range
#strand: + or -
#metadata: anything

#require(biomaRt)
# set directories #
# parent directory
home_dir <- getwd()

# 2ndStruc save dir
save_dir <- paste0(home_dir,"/../Chapter_4/call_sites")
dir.exists(save_dir)
# Filter table
# Download from https://github.com/ZW-xjtlu/rBS2ndStructure/blob/master/data/Struc_mm10.rda
load(paste0(save_dir,"/../Struc_mm10.rda"))

## --- read in raw data files --- ##
list_of_anno_files <-c('SRR8170379')#, 'SRR8170378','MF_rep2') #c('G1', 'G2', 'G3', 'G4', 'MF_rep1', 'MF_rep2',
                       #'pMF_rep1', 'pMF_rep2', 'SRR8170377', 'SRR8170378', 'SRR8170379', 'SRR8170380')

get_2ndfilter_individualFiles <- function(test_m5C) {
  # Import, format, merge dataset to mm10 reference
  name <- test_m5C
  rBS_df_OG <- read.csv(paste0(save_dir, "/", test_m5C,'_Genome10xCall_signalFilter.txt'), sep = "\t")
  rBS_df_OG$site_loc <- paste0(rBS_df_OG$X.SeqID, "_", rBS_df_OG$refPos)
  
  # create column for absolute position reference, (subtract on neg strans, add on pos strand)
  # format chrom to match ref file
  rBS_gr <- GRanges(seqnames = paste0("chr",gsub("MT","M",rBS_df_OG$X.SeqID)),
                    strand = rBS_df_OG$refStrand,
                    ranges = IRanges(start = rBS_df_OG$refPos, width = 1), 
                    site_loc = rBS_df_OG$site_loc)
  # Add filter controls if needed 
  #pos_cont <- tail(Struc_mm10, n=1)
  #pos_cont <- data.frame(pos_cont)
  #pos_cont <- makeGRangesFromDataFrame(pos_cont)
  #rBS_final <- c(rBS_gr, pos_cont)
  
  # perform filter
  cat(paste(test_m5C,": Before filter",nrow(data.frame(rBS_df_OG)),"\n"))
  rBS_gr_retained <- rBS_gr[!rBS_gr %over% Struc_mm10]
  cat(paste(test_m5C,": After filter",length(seqnames(rBS_gr_retained)),"\n"))
  # write sites removed and retained files (subset of og file)
  #removed <- rBS_df_OG[rBS_df_OG$group %in% rBS_gr_retained$group,]
  retained <- rBS_df_OG[rBS_df_OG$site_loc %in% rBS_gr_retained$site_loc,]
  write.table(retained, paste0(save_dir,"/", name, "_strucFilter.txt"),
              quote=FALSE, sep="\t", row.names=FALSE)

  
}
#test
#get_2ndfilter_individualFiles(test_m5C)
#apply to file list
lapply(list_of_anno_files,get_2ndfilter_individualFiles)



