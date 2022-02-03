#!/bin/bash

module reset
module load UMI-tools

cd $2/result/$1

umi_tools dedup --paired --stdin=$1_meRanGh_genomeMap_sorted.bam --stdout=$1_meRanGh_genomeMapTrim_dedup.bam --output-stats=$1_meRanGh_genomeMap_dedupStats.txt
umi_tools group --paired --output-bam --stdin=$1_meRanGh_genomeMap_sorted.bam --stdout=$1_meRanGh_genomeMapTrim_dedupGroup.bam --group-out=$1_meRanGh_genomeMap_UMIgrouped.tsv

