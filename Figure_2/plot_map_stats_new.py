import os, re
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# %%
# methods
"""
Code to plot mapping efficiencies from SAMtools output files.
"""

# Read C_counts from fasta reference
def countReads_from_flagstat(input_flagstat):
    with open(input_flagstat, "r") as input:
        input_stats = input.readline()
        input_stats = int(input_stats.split()[0])

    return input_stats

def countReads_from_fqCount(input_fqCount):
    with open(input_fqCount, "r") as input:
        fqCount = int(input.readline())

    return fqCount

# %%

"""
Code will search for files with prefix "file" and match statistics files 
Args:
    meRan*count.txt: SAMtools flagstat output
    fq.gz_count.txt: A single line text file only containing read count (zcat $file | wc -l)
Return:
    Dataframe of all samples in glob list ready for easy plotting
"""

plot_dict = {}
norm_compare = {}
raw_compare = {}
for file in glob.glob("**/*_meRan*_count.txt", recursive=True):  # remove '**/' if files are local to cwd
    name = os.path.basename(file).split("_count.txt")[0]
    baseName = re.sub(r'_meRan.*$','', name)
    fqCount_name = baseName + '_1_Window_4bpTrim.fq.gz_count.txt'

    try:
        fqCount_file = glob.glob("**/" + fqCount_name, recursive=True)
        print(file)
        plot_dict[name] = [countReads_from_flagstat(file) / (countReads_from_fqCount(fqCount_file[0])*2)]
    except:
        print("FQ file count file not found")
tall_df = pd.DataFrame.from_dict(plot_dict).T

# %%
def format_df(df):
    df['name'] = df.index.astype(str).str.replace('_meRan.*$','')
    df['map_approach'] = np.where(df.index.str.contains("_GenomeMulti"), "GN2TX",
                        np.where(df.index.str.contains("_genomeMap"), "Genome",
                        np.where(df.index.str.contains("_TXMap"), "Trans", "Not found")))
    #df['map_approach'] = np.where(df['map_approach'] == True, "Trans", "Genome")
    df['group'] = df['name'].str.replace(r'817.*$|_.*$', '')
    df = df.sort_values('map_approach')
    return df
test_df = format_df(tall_df)
#%%
# Add multimapped-mapped to true map
test = pd.pivot(test_df.reset_index(), index='name', columns='map_approach', values=0)
test['GN_plus_GNmulti2TX'] = test.Genome + test.GN2TX
test = test.drop(['GN2TX','Not found'], axis=1)

# reshape to long format
plot_df = pd.melt(test.reset_index(), id_vars='name')
plot_df['group'] = plot_df['name'].str.replace(r'817.*$|_.*$', '')
#%%
plot_df.to_csv("map_stats.csv")
# %%

# All samples' mapping rate
order_x = ['G5', 'G1', 'G2', 'G3', 'G4', 'MFRNAseq', 'MF', 'pMFRNAseq', 'pMF', 'SRR']
#order_hue = ['Genome','Trans','GN2TX','TX2GN']
sns.set(rc={'figure.figsize':(10,5)})
sns.despine()
sns.set_style("ticks")

p = sns.barplot(data=plot_df, x='group', y='value', hue='map_approach', edgecolor="black", linewidth=2,
                palette='bright', order=order_x)

plt.xticks(rotation=90)
p.set(yticks=np.arange(0, 1.1, 0.5))

plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
plt.savefig('Perc_mapped_reads.png', bbox_inches='tight', dpi=400, transparent=True)

plt.show()
