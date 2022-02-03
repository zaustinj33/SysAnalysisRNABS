#!/bin/bash

## Input txt file should be the desired output file prefix
## ex if input files are MF_batch1_rep1 MF_batch2_rep1 MF_batch1_rep2 MF_batch2_rep2 ...
## list.txt should be MF_rep1 MF_rep2 ....
## manually delete 3rd replicate in samtools merge if it doesn't exist 

NAME=$(echo $1 | cut -d'_' -f 1) # get base name of split files
REP=$(echo $1 | cut -d'_' -f 2) # get replicate name of split files

cd $2
echo '$NAME $REP'

module reset
module load SAMtools
#1 merge bam files
samtools merge -@ 40 -f -r $2/result/$1/$1_Merged_final.bam $2"/result/"$NAME"_batch1_"$REP"/"$NAME"_batch1_"$REP"_meRanGh_genomeMap_sorted.bam" $2"/result/"$NAME"_batch2_"$REP"/"$NAME"_batch2_"$REP"_meRanGh_genomeMap_sorted.bam" $2"/result/"$NAME"_batch3_"$REP"/"$NAME"_batch3_"$REP"_meRanGh_genomeMap_sorted.bam"

#2 sort merged file
samtools sort -@ 40 $2/result/$1/$1_Merged_final.bam > $2/result/$1/$1_meRanGh_genomeMapMerged_sorted.bam
rm $2/result/$1/$1_Merged_final.bam

#3 index sort-merged file
samtools index -@ 40 $2/result/$1/$1_meRanGh_genomeMapMerged_sorted.bam
