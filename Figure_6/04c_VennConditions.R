library(grid)
library(VennDiagram)
library(scales)
library(gplots)

## Venn diagram to compare m5C replicates and overlap

# set directories #
home_dir <- paste0(getwd(), '/Chapter_5')
setwd(home_dir)

union_files <- list.files(pattern = "allOverlapDepth.csv", full.names = T)  # Union
union_files[1] <- "./G_allUnion.csv"

generate_gene_list <- function(file){
  overlap <- read.csv(file)
  name <-  gsub("all.*$","",basename(file))

  gene_list <- overlap$group
  return(gene_list)
}
input <- lapply(union_files[-4], generate_gene_list)

pdf(file=paste0("G_MF_pMF_venn.pdf"))
venn_table <- venn(input, intersections = TRUE)
dev.off()

intersections <- attr(venn_table,"intersections")
union <- unlist(intersections)
overlap <- intersections[[1]]

