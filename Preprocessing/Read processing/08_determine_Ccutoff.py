import glob
import pandas as pd
import os, sys, re
import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns

#%%

#1) filter df by basic requirements, and write output to file
def filter_call_file(input_file):
    name = re.sub("_Genome10xCall*_annotate.txt", "", input_file)  # sample name
    long_name = re.sub("_annotate.txt", "", input_file)  # write file
    sys.stdout.write("\n")
    sys.stdout.write("Processing " + os.path.basename(name) + '\n')
    sys.stdout.write("\n")
    df = pd.read_csv(input_file, sep='\t', low_memory=False)  # file

    # Filters:
    # 1) coverage >= 20
    # 2) C count >= 3
    # 3) methylation rate (count/coverage) >= 0.1
    # 4) reference base = C
    # 5) must be assigned to a feature
    final_df = df[((df['cov'] >= 20) & (df['methRate'] >= 0.1) & (df['C_count'] >= 3) &
                   (df['state'] == 'M'))]
    
    # Cutoff for C_count > 3 only
    #final_df = df[((df['methRate'] >= 0.1) & (df['C_count'] >= 3) & (df['state'] == 'M'))]

    final_df.to_csv(long_name+"_basicFilter.txt", index=False, sep='\t')

    filtered_call = final_df[(df['gene'] != 'gene:no_feature;')]

    return filtered_call, name


# %%

def gini(array):
    """Calculate the Gini coefficient of a numpy array.
    usage: for a list of genes sorted by number of m5C sites, calculate Gini:
        sum(2i - n - 1) * N_i / (n * sum(N_i))
        i = index of gene
        N_i = number of m5C sites in the ith gene
        n = total number of unique genes
    """
    """
    # Method 1) Naive Gini
    # based on:
    # http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.html
    # oliviaguest@github.com

    # All values are treated equally, arrays must be 1d:
    array = array.flatten()

    # Values must be sorted:
    array = np.sort(array)

    # Index per array element:
    index = np.arange(1, array.shape[0] + 1)

    # Number of array elements:
    n = array.shape[0]

    # Gini coefficient:
    return (np.sum((2 * index - n - 1) * array)) / (n * np.sum(array))
    """

    # Method 2) adapted from Huang et al
    sorted_list = sorted(array)
    height, area = 0, 0
    for value in sorted_list:
        height += value
        area += height - value / 2.
    fair_area = height * len(array) / 2.
    return (fair_area - area) / fair_area


#%%
Ccutoff = [None, 15, 8, 5, 4, 3, 2, 1]
def calc_Gini_from_bam(filtered_call, name, cutoff):

        sys.stdout.write("Testing C-cutoff of: " + str(cutoff))
        df = filtered_call
        df_gene_counts = df.gene.value_counts()
        gini_value = round(gini(np.asarray(df_gene_counts, dtype=np.float64)),6)

        sys.stdout.write(" GINI: " + str(gini_value) + '\n' + "Total counts: " + str(sum(df_gene_counts))+'\n')
        return gini_value


# %%
cutoff_dict = {}
sample_list = ['MF_rep2'] #['SRR8170377', 'SRR8170378', 'SRR8170379', 'SRR8170380', 'G1','G2','G3','G4', 'MF_rep1', 'MF_rep2','pMF_rep1', 'pMF_rep2']
#['Gall', 'G1','G1_1A', 'G1_2A','G2', 'G2_1B', 'G2_2B','G3','G3_1C', 'G3_2C', 'G4', 'G4_1D', 'G4_2D']
for sample in sample_list:
    cutoff_dict[sample] = {}
    for cutoff in Ccutoff:
        gini_value = None
        if cutoff == None:  # no c-cutoff
            target_file = sample+"_Genome10xCall_annotate.txt"
        else:
            target_file = sample+"_Genome10xCall_"+str(cutoff)+"_Cutoff_annotate.txt"
        print(target_file)
        for file in glob.glob("**/"+target_file, recursive=True):
            print("found")
            try:
                filtered_df, input_name = filter_call_file(file)
                gini_value = calc_Gini_from_bam(filtered_df, input_name, cutoff)
                cutoff_dict[sample][cutoff] = gini_value


                gene_counts = filtered_df.gene.value_counts()

                ## Set up histplot
                num_one = len(gene_counts[gene_counts == 1])
                if max(gene_counts) > 25:
                    xmax = max(gene_counts)
                else:
                    xmax = 25

                sns.histplot(data=gene_counts, color='green', bins=range(1,max(gene_counts)))#,log=True)
                plt.xlim(0,xmax)
                plt.yscale('log')
                #plt.text(xmax-xmax*0.1, num_one-num_one*0.1,"Gini: " + str(round(gini_value,4)) +
                 #        "\nC-cutoff: " + str(cutoff), va='center', ha='right')
                plt.xlabel("Sites per gene")
                plt.ylabel("Count")
                plt.title(sample)
                #
                #plt.savefig(os.getcwd() + "/" + sample + "_cutoff_" + str(cutoff)
                 #           + "_sitesPerGene_histGini.png", bbox_inches='tight', dpi=400, transparent=True)
                #plt.show()
                plt.close()

            except ValueError:
                print("no sites found in file \n")
                cutoff_dict[sample][cutoff] = gini_value
                pass
            except TypeError:
                print("weird value in column?\n")
                cutoff_dict[sample][cutoff] = gini_value
                pass




#%%

with open(os.getcwd()+'/Gall_Ginivalues.json', 'w') as convert_file:
    convert_file.write(json.dumps(cutoff_dict))
