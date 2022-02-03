#!/bin/bash

GNconvRate=$(cat $2/Code/$1/Ch3_code/$1_GN_conv_rate.txt)

module reset
module load SAMtools

## Genome
# Call sites
meRanCall -p 40 -o $2/CallResult/$1/$1_Genome10xCall.txt -bam $2/result/$1/$1_meRanGh_genomeMapMerged_sorted.bam -f /projects/epigenomicslab/Annotation/mm10.for.RNABS.fa -mBQ 30 -gref -rl 150 -sc 10 -cr $GNconvRate -fdr 0.01 -mr 0.00001 -mcov 10

# Annotate sites
sed -i 's/chrM/chrMT/g' $2/CallResult/$1/$1_Genome10xCall.txt
sed -i 's/chr//' $2/CallResult/$1/$1_Genome10xCall.txt 

meRanAnnotate -t $2/CallResult/$1/$1_Genome10xCall.txt -ensGTF -g /projects/epigenomicslab/Annotation/Mus_musculus.GRCm38.96.gtf -o $2/CallResult/$1/$1_Genome10xCall_annotate.txt -f 'gene'
