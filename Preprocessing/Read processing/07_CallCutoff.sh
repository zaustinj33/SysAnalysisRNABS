#!/bin/bash

module reset
module load SAMtools

mkdir -p \$TMPDIR/ZJ_tmp

## Genome
# Call sites
# Note: very I/O intensive. Try to create files on job node, then move to local

for i in 1 2 3 4 5 8 15; do
    echo "cutoff: \$i"
    samtools index $2/result/$1/$1_meRanGh_genomeMapMerged_\"\$i\"_Ccutoff_PE.bam
    meRanCall -p 40 -o \$TMPDIR/ZJ_tmp/$1_Genome10xCall_\"\$i\"_Cutoff.txt -bam $2/result/$1/$1_meRanGh_genomeMapMerged_\"\$i\"_Ccutoff_PE.bam -f /projects/epigenomicslab/Annotation/mm10.for.RNABS.fa -mBQ 30 -gref -rl 150 -sc 10 -cr 1 -fdr 0.01 -mr 0.00001 -mcov 10
    mv \$TMPDIR/ZJ_tmp/$1_Genome10xCall_\"\$i\"_Cutoff.txt $2/CallResult/$1/$1_Genome10xCall_\"\$i\"_Cutoff.txt

    sed -i 's/chrM/chrMT/g' $2/CallResult/$1/$1_Genome10xCall_\"\$i\"_Cutoff.txt
    sed -i 's/chr//' $2/CallResult/$1/$1_Genome10xCall_\"\$i\"_Cutoff.txt 

    meRanAnnotate -t $2/CallResult/$1/$1_Genome10xCall_\"\$i\"_Cutoff.txt -ensGTF -g /projects/epigenomicslab/Annotation/Mus_musculus.GRCm38.96.gtf -o $2/CallResult/$1/$1_Genome10xCall_\"\$i\"_Cutoff_annotate.txt -f 'gene'
done

# No cutoff
samtools index $2/result/$1/$1_meRanGh_genomeMapMerged_sorted.bam
meRanCall -p 40 -o \$TMPDIR/ZJ_tmp/$1_Genome10xCall.txt -bam $2/result/$1/$1_meRanGh_genomeMapMerged_sorted.bam -f /projects/epigenomicslab/Annotation/mm10.for.RNABS.fa -mBQ 30 -gref -rl 150 -sc 10 -cr 1 -fdr 0.01 -mr 0.00001 -mcov 10
mv \$TMPDIR/ZJ_tmp/$1_Genome10xCall.txt $2/CallResult/$1/$1_Genome10xCall.txt


# Annotate sites
sed -i 's/chrM/chrMT/g' $2/CallResult/$1/$1_Genome10xCall.txt
sed -i 's/chr//' $2/CallResult/$1/$1_Genome10xCall.txt

meRanAnnotate -t $2/CallResult/$1/$1_Genome10xCall.txt -ensGTF -g /projects/epigenomicslab/Annotation/Mus_musculus.GRCm38.96.gtf -o $2/CallResult/$1/$1_Genome10xCall_annotate.txt -f 'gene'

" > $1_meRanCallCutoff.sbatch

sbatch $1_meRanCallCutoff.sbatch
