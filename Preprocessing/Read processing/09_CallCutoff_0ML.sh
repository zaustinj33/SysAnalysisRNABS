#!/bin/bash

GNconvrate=$(cat $2/Code/$1/Ch3_code/$1_GN_conv_rate.txt)
mkdir -p \$TMPDIR/ZJ_tmp

# Usage: bash 07_... [name of sample] [working directory] [cutoff]
# $2/All_m5C.bed is output from 02b_Create_raw_union_file.R

## Genome
meRanCall -p 40 -o \$TMPDIR/ZJ_tmp/$1_Genome10xCall_$3_Cutoff_0ML.txt -bam $2/result/$1/$1_meRanGh_genomeMapMerged_$3_Ccutoff_PE.bam -f /projects/epigenomicslab/Annotation/mm10.for.RNABS.fa -mBQ 30 -gref -rl 150 -sc 10 -cr 1 -mr 0 -mcov 10 -regions $2/All_m5C.bed
mv \$TMPDIR/ZJ_tmp/$1_Genome10xCall_$3_Cutoff_0ML.txt $2/CallResult/$1/$1_Genome10xCall_$3_Cutoff_0ML.txt
