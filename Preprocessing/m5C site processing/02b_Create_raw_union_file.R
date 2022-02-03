library(purrr)
library(dplyr)
library(GenomicRanges)

setwd(paste0(getwd(),"/Chapter_5"))

# Find files to merge

# change suffix to find appropriate files
suffix <- "_0ML.txt"  #0ML depth sites
suffix <- "_strucFilter.txt"  # sites after secondary structure filter
path_fdr <- "/../Chapter_4/call_sites/"
#path_struc <- "/results/"

mRNAfiles <- list.files(path=paste0(getwd(),path_fdr), pattern = suffix,full.names = T, recursive = T)

mRNA_dfs <- lapply(mRNAfiles, function(x) {
  # read in raw files and grab only relevant columns to merge together
  y <- read.csv(x, comment.char = '',header = T,sep='\t')[,c(1,2,5,6,7,3,18)] #21
  name <- gsub(suffix,"",basename(x))

  # add name from column for easy merging
  names(y)[-6:-8] <- paste0(names(y)[-6:-8],"_",name)
  
  
  #20x, 3m5C filter for replicates
  # comment out if preparing depth file 
  #y <- y %>%
  #  filter_if(grepl("C_count",colnames(y)), any_vars(. >= 3)) %>%
  #  filter_if(grepl("cov",colnames(y)), any_vars(. >= 20))


  cat(name,": ",nrow(y),"\n")
  y$group <- paste(y$X.SeqID, y$refPos)
  y[,-1:-2]
})


# Combined sites dataframe
merge.total <- Reduce(function(...) merge(..., by = c('group','seqContext','refStrand'), all = T), mRNA_dfs) #'gene'
merge.total <- merge.total[!duplicated(merge.total$group),]
merge.total$group <- gsub("chr","",merge.total$group)
merge.total$group <- gsub("M","MT",merge.total$group)
merge.total$chrom <- gsub(" .*","",merge.total$group)
merge.total$refPos <- as.integer(gsub(".* ","",merge.total$group))
#merge.total <- merge.total %>%
 # filter_if(grepl("C_count",colnames(merge.total)), any_vars(. >= 3)) %>%
  #filter_if(grepl("cov",colnames(merge.total)), any_vars(. >= 20))


#merge.total.final <- merge.total
#write.csv(merge.total, "allm5C_libraries.csv",row.names = F)
write.csv(merge.total, "allm5C_libraries_filtered.csv",row.names = F)

## Write bed file for calling sites at depth
all_bed <- merge.total %>% select(chrom, refPos)
all_bed$chrom <- paste0("chr", all_bed$chrom)
all_bed$chrom <- gsub("chrMT", "chrM", all_bed$chrom)
all_bed$end <- all_bed$refPos+1
all_bed$refPos <- all_bed$refPos-1
write.table(all_bed, "All_m5C.bed", quote = F, sep='\t', row.names = F, col.names = F)

## Filter unfiltered df by filtered df to get depth at all sites in the

unfiltered <- read.csv("allm5C_libraries.csv", comment.char = '',header = T)
filtered <- read.csv("allm5C_libraries_filtered.csv", comment.char = '',header = T)

filtered_depth <- subset(merge.total,group %in% filtered$group)
filtered_depthAnno <- left_join(filtered_depth, filtered, by="group")
write.csv(filtered_depth, "allm5C_libraries_filteredDepth.csv",row.names = F)
