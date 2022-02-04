#!/bin/bash

# Adapter trim #

module reset
module load fastp

cd $2/raw_data/$1
# 6bp trim
fastp -w 16 -f 6 -t 6 -i $1_1_qual.fq.gz -I $1_2_qual.fq.gz --out1 $1_1_Window.fq.gz --out2 $1_2_Window.fq.gz --failed_out $1_failed_Window.fq.gz -j $1_failed_Window.json -h $1_failed_Window.html

module reset
module load HISAT2
gunzip $1_1_Window.fq.gz
gunzip $1_2_Window.fq.gz

echo 'Start meRanGh'
meRanGh align -o $2/result/$1/ -un -ud $2/raw_data/$1/meRanGhUnaligned -f $2/raw_data/$1/$1_2_Window.fq -r $2/raw_data/$1/$1_1_Window.fq -t 40 -fmo -S $2/result_new/$1/$1_meRanGh_genomeMap.sam -ds -MM -id /projects/epigenomicslab/Annotation/mm10_meRanGh/BSgenomeIDX -GTF /projects/epigenomicslab/Annotation/mm10.ensGene.for.RNABS.gtf

echo 'Start meRanT'
meRanT align -o $2/result/$1/ -un -ud $2/raw_data/$1/meRanTUnaligned -f $2/raw_data/$1/$1_2_Window.fq -r $2/raw_data/$1/$1_1_Window.fq -t 40 -fmo -S $2/result_new/$1/$1_meRanT_TXMap.sam -ds -ra -MM -i2g /projects/epigenomicslab/Annotation/mouse.rna.ensemble.ERCC.map -x /projects/epigenomicslab/Annotation/mm10_meRanT/mm10.for.RNABS_C2T
echo 'Finished meRanGh'

rm $1_1_Window.fq
rm $1_2_Window.fq

cd $2/raw_data/$1/meRanTUnaligned
gunzip $1_2_Window_unmapped.fq.gz
gunzip $1_1_Window_unmapped.fq.gz

## Map multimapeprs
echo 'Start meRanGh'
meRanGh align -fmo -o $2/result/$1/ -f $2/raw_data/$1/meRanTUnaligned/$1_2_Window_unmapped.fq -r $2/raw_data/$1/meRanTUnaligned/$1_1_Window_unmapped.fq -t 40 -S $2/result_new/$1/$1_meRanGh_TransMulti_to_GN.sam -ds -MM -id /projects/epigenomicslab/Annotation/mm10_meRanGh/BSgenomeIDX -GTF /projects/epigenomicslab/Annotation/mm10.ensGene.for.RNABS.gtf
echo 'Finished meRanGh'

rm $1_2_Window_unmapped.fq
rm $1_1_Window_unmapped.fq

cd $2/raw_data/$1/meRanGhUnaligned
gunzip $1_2_Window_unmapped.fq.gz
gunzip $1_1_Window_unmapped.fq.gz

echo 'Start meRanT'
meRanT align -fmo -o $2/result/$1/ -f $2/raw_data/$1/meRanGhUnaligned/$1_2_Window_unmapped.fq -r $2/raw_data/$1/meRanGhUnaligned/$1_1_Window_unmapped.fq -t 40 -S $2/result_new/$1/$1_meRanT_GenomeMulti_to_TX.sam -ds -ra -MM -i2g /projects/epigenomicslab/Annotation/mouse.rna.ensemble.ERCC.map -x /projects/epigenomicslab/Annotation/mm10_meRanT/mm10.for.RNABS_C2T

rm $1_2_Window_unmapped.fq
rm $1_1_Window_unmapped.fq