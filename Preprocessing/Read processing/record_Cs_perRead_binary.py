import glob, os, json
import pysam

## Note: pysam does not work in Windows!

#%%

def countCs_from_BAM(bam_input):
    #print(bam_input)
    C_count_dict = {'no C':0,'with C':0}
    C_count_list = []
    bamfile = pysam.AlignmentFile(bam_input,"rb", check_sq=False)

    # check each read for C content, write to new file if passed filter
    for read in bamfile.fetch():
        sequence = read.query_sequence.upper()
        C_count = 'null'
        try:
            if read.get_tag('YR') == 'G2A':  # read pair copy
                pass
            elif read.get_tag('YG') == 'C2T':  # Aligned to positive strand
                C_count = float(sequence.count('C'))
            elif read.get_tag('YG') == 'G2A':  # Aligned to negative strand
                C_count = float(sequence.count('G'))
            else:
                print("tag not found")
                pass
        except KeyError:
            pass

        #C_count_list.append(C_count)
        # bin C_count per read
        if C_count == 'null':
            pass
        elif C_count == 0:
            C_count_dict['no C'] += 1
        elif C_count >= 1:
            C_count_dict['with C'] += 1
        else:
            pass

    return C_count_dict

#%%
C_count_dict = {}
df_list = []
name_list = []

# Important files to count: G* (MT), *Merged (pMF, MF), SRR* (Huang), *RNAseq* (pMF and MF RNAseq)
important_ext = ["G*genomeMap_dedup.bam"]
#["G*genomeMap_sorted.bam", "*genomeMapMerged_sorted.bam", "SRR*genomeMap_sorted.bam","*RNAseq*genomeMap_sorted.bam"]
important_files = []

for ext in important_ext:
    for file in glob.glob("**/"+ext, recursive=True):
        important_files.append(file)

print(important_files)
for file in important_files:
    print(file)
    name = os.path.basename(file).split("*genomeMap_dedup.bam")[0]
    name_list.append(name)
    hist_list = countCs_from_BAM(file)
    df_list.append(hist_list)

for i in range(len(name_list)):
    C_count_dict[name_list[i]] = df_list[i]

with open(os.getcwd()+'/C_Count_hist_binary_Dedup.json', 'w') as convert_file:
    convert_file.write(json.dumps(C_count_dict))


