#!/bin/bash

mkdir -p $2/Code/$1/TecPap
mkdir -p $2/mapped_files/$1
mkdir -p $2/raw_data/$1
mkdir -p $2/mapped_files/$1/WithOverlap
mkdir -p $2/mapped_files/$1/NoOverlap

pushd $2/Code/$1/TecPap
RAW1=$(find $2 -name $1_1.fq.gz)
RAW2=$(find $2 -name $1_2.fq.gz)
UMI=$(find $2 -name $1_umi.fq.gz)

# Quality analysis breakdown for RNA-BS reads, outputs reads failed at each step to seperate files for downstream analysis #


cd $2/raw_data/$1

#Add UMI info to headers
module load UMI-tools
umi_tools extract --bc-pattern=NNNNNNNN --stdin=$UMI --read2-in=$RAW1 --stdout=$1_1_umi.fq.gz --read2-stdout
umi_tools extract --bc-pattern=NNNNNNNN --stdin=$UMI --read2-in=$RAW2 --stdout=$1_2_umi.fq.gz --read2-stdout

echo 'Starting quality analysis'

## Step by step quality analysis to determine source of bad reads
module reset
module load fastp

# Too short or no mate, remove adapters, remove polyx tails
fastp -w 16 -Q -l 50 --trim_poly_x --poly_x_min_len 10 -i $1_1.fq.gz -I $1_2.fq.gz --out1 $1_1_trim.fq.gz --out2 $1_2_trim.fq.gz --failed_out $1_failed_length.fq.gz -j $1_failed_length.json -h $1_failed_length.html --detect_adapter_for_pe --overlap_diff_percent_limit 25

# Remove low quality bases
fastp -w 16 -q 25 -i $1_1_trim.fq.gz -I $1_2_trim.fq.gz --out1 $1_1_qual.fq.gz --out2 $1_2_qual.fq.gz --failed_out $1_failed_qual.fq.gz -j $1_failed_qual.json -h $1_failed_qual.html --detect_adapter_for_pe

# re-attach UMIs to end of read with '+' separator
zcat $1_1_qual.fq.gz | sed -E 's/([A-Z]{8,} )(.*)/\2+\1/' | gzip > $1_1_qualFormatted.fq.gz
zcat $1_2_qual.fq.gz | sed -E 's/([A-Z]{8,} )(.*)/\2+\1/' | gzip > $1_2_qualFormatted.fq.gz

## get mbias data
# map files
module reset
module load Bismark
bismark -p 4 --genome /projects/epigenomicslab/Annotation/mm10_bismark/ -1 $1_2_qual.fq.gz -2 $1_1_qual.fq.gz -o $2/mapped_files/$1
cd $2/mapped_files/$1

module load SAMtools
samtools sort -@ 8 $1_1_qualFormatted_bismark_bt2_pe.bam > $1_sort.bam
samtools index $1_sort.bam

module reset
module load UMI-tools
umi_tools dedup --paired --umi-separator=+ -I $1_sort.bam --output-stats=$1_deduplicated -S $1_deduplicated.bam
umi_tools group --paired --umi-separator=+ -I $1_sort.bam --group-out=$1_groups.tsv --output-bam -S $1_groups.bam

rm $1_sort.bam

samtools sort -n -@ 8 $1_deduplicated.bam > $1_deduplicated_sort.bam
samtools index $1_deduplicated_sort.bam

module reset
module load Bismark
bismark_methylation_extractor -p --parallel 10 --gzip --no_overlap --genome_folder /projects/epigenomicslab/Annotation/mm10_bismark/ $1_2_Window_bismark_bt2_pe.bam -o $2/mapped_files/$1/NoOverlap

rm $1_deduplicated.bam
