library(grid)
library(VennDiagram)
library(scales)
library(gplots)

## Venn diagram to compare m5C replicates and overlap

# set directories #
home_dir <- paste0(getwd(), '/Chapter_5')
setwd(home_dir)

union_files <- list.files(pattern = "allOverlapRepDepth.csv", full.names = T)  # Union
#union_files[1] <- "./G_allUnion.csv"

generate_gene_list <- function(file){
  overlap <- read.csv(file)
  name <-  gsub("all.*$","",basename(file))

  gene_list <- overlap$group
  return(gene_list)
}
input <- lapply(union_files, generate_gene_list)

pdf(file=paste0("G_SRR_venn.pdf"))
venn_table <- venn(input, intersections = TRUE)
dev.off()

intersections <- attr(venn_table,"intersections")
union <- unlist(intersections)

length(intersections[[1]]) <- length(intersections[[2]])
length(intersections[[3]]) <- length(intersections[[2]])
GO_input <- as.data.frame(intersections, check.rows = FALSE, col.names = c("MT_only", "Huang_only", "Overlap"))

write.table(GO_input,"GOinput_huang_MT_overlap.txt", sep='\t', quote=F, row.names=F)

## comparing venn intersections
all_m5C <- read.csv("allm5C_libraries_filteredDepthAnno.csv")
all_highConf_sites <- subset(all_m5C, group %in% union)
write.csv(all_highConf_sites, "allm5C_highConf_sites.csv", row.names = F)
MT_only <- subset(all_highConf_sites, group %in% intersections[[1]])
Huang_only <- subset(all_highConf_sites, group %in% intersections[[2]])
MT_Huang_ol <- subset(all_highConf_sites, group %in% intersections[[3]])


all_m5C_filtered <- all_highConf_sites %>%
  filter_if(grepl("C_count",colnames(all_m5C)), any_vars(. >= 3)) %>%
  filter_if(grepl("cov",colnames(all_m5C)), any_vars(. >= 20))

Huang_sites <- read.csv("Huang_m5C_muscle_published.csv")
Huang_sites$group <- paste(Huang_sites[,1], Huang_sites$Position)

pdf(file=paste0("Huang_SRR_venn.pdf"))
venn_table <- venn(list(Huang_sites$gene_name[Huang_sites$if.m5C.site == 'TRUE'], Huang_only$gene_name), intersections = TRUE)
dev.off()
