#%%
import pysam
import glob, re, sys

#%%

Ccutoff_list = {1:0, 2:0, 3:0, 4:0, 5:0, 8:0, 15:0}

# Define functions
"""Count occurances of C's in bam file, remove them, write clean bam file to new file. """
def count_C(samfile):
    name = re.sub("_sorted.bam",'',samfile)
    #name = re.sub("_Map_sorted.bam",'',samfile)
    #name = re.sub("^subsample_",'',name)
    sys.stdout.write("processing:\t" + samfile)
    removed_read_count = 0
    bamfile = pysam.AlignmentFile(samfile,"rb", check_sq=False, threads=4)

    # Open cutoff files for writing
    for k in Ccutoff_list.keys():
        Ccutoff_list[k] = pysam.AlignmentFile(name+"_"+str(k)+"_Ccutoff_PE.bam", "wb", template=bamfile)
    # check each read for C content, write to new file if passed filter
    for read in bamfile.fetch():
        C_count = 'null'
        sequence = read.query_sequence.upper()
        try:
            #if read.get_tag('YR') == 'G2A':  # read pair copy
             #   pass
            if read.get_tag('YG') == 'C2T':  # Aligned to positive strand
                C_count = float(sequence.count('C'))
            elif read.get_tag('YG') == 'G2A':  # Aligned to negative strand
                C_count = float(sequence.count('G'))
            else:
                print("tag not found")
                pass
        except KeyError:
            print("tag not found")
            break

        if C_count == 'null':
            pass
        elif C_count <= max(Ccutoff_list):
            for k,v in Ccutoff_list.items():
                if C_count <= k:
                    v.write(read)
        else:
            removed_read_count += 1
    bamfile.close()
    for k,v in Ccutoff_list.items():
        v.close()
    sys.stdout.write(str(removed_read_count) + " reads removed \n")
    return removed_read_count


#%%
#for bam_file in glob.glob("**/G*_genomeMapMerged_sorted.bam", recursive=True):
 #   count_C(bam_file)
 
#%%
#for bam_file in glob.glob("**/MF_rep2*_genomeMapMerged_sorted.bam", recursive=True):
 #   count_C(bam_file)

#%%
for bam_file in glob.glob("**/SRR*80*_genomeMap_sorted.bam", recursive=True):
    count_C(bam_file)

#%%
#for bam_file in glob.glob("**/*_genomeMap_dedup.bam", recursive=True):
 #   count_C(bam_file)

#%%
#for bam_file in glob.glob("**/*RNAseq*_genomeMap_sorted.bam", recursive=True):
    #count_C(bam_file)

