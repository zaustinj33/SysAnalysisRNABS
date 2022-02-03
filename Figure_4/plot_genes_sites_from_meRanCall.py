import glob
import pandas as pd
import os, sys, re
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

#%%
headers = ['#SeqID', 'refPos', 'strand', 'Base', 'cov', 'C_count','methRate' ]
methRate_counts = {}
methRate_counts_filtered = {}
gene_counts = {}
sites_per_gene = {}
for file in glob.glob("**/*_Genome10xCall_annotate.txt", recursive=True):  #*_Genome10xCall*_annotate.txt", recursive=True):
    name = re.sub("_Genome10xCall_annotate.txt", "", os.path.basename(file))
    print(file)
    df = pd.read_csv(file, sep='\t', low_memory=False)  # file
    #df = df.iloc[:,:7]
    #df.columns = headers

    #1) basic filter
    df_filter = df[((df['cov'] >= 20) & (df['methRate'] >= 0.1) & (df['C_count'] >= 3))]
    sites_per_gene_sample = df.gene.value_counts()
    #2) genes per library
    gene_counts[name] = len(sites_per_gene_sample)
    #3) sites per gene
    sites_per_gene[name] = sites_per_gene_sample
    #4) methylation rate
    methRate_counts[name] = df['methRate']
    methRate_counts_filtered[name] = df_filter['methRate']

#%%
order_x = ['G1', 'G2', 'G3', 'G4', 'MF', 'pMF', 'SRR']

# plot genes per sample
gene_counts_df = pd.DataFrame(gene_counts.items(), columns=['name', 'counts'])
gene_counts_df['group'] = gene_counts_df['name'].str.replace(r'817.*$|_.*$', '')  # add levels

sns.barplot(data=gene_counts_df, x='group', y='counts', palette='bright', order=order_x)
sns.swarmplot(data=gene_counts_df, x='group', y='counts', color='black', order=order_x)

plt.xticks(rotation=90)
plt.savefig("genes_per_lib.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()
plt.close()
#%%

# plot sites per gene
sites_df = pd.DataFrame(sites_per_gene)
sites_melt = sites_df.melt()
#%%
index_order = ['G1', 'G2','G3', 'G4', 'MF_rep1',
               'MF_rep2', 'pMF_rep1',
               'pMF_rep2', 'SRR8170377','SRR8170378', 'SRR8170379', 'SRR8170380']

sns.violinplot(data=sites_melt, x='variable', y='value', color='black', order=index_order,
               cut=0)  # inner='point'
plt.xticks(rotation=90)
plt.yscale('log')
plt.savefig("sites_per_gene.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()
plt.close()

#%%
#methRate_df
methRate_df = pd.DataFrame(methRate_counts).melt()
methRate_df['group'] = methRate_df['variable'].str.replace(r'817.*$|_.*$', '')
#%%
methRate_df_filtered = pd.DataFrame(methRate_counts_filtered).melt().dropna()
methRate_df_filtered['group'] = methRate_df_filtered['variable'].str.replace(r'817.*$|_.*$', '')

#%%
palette_bright = sns.color_palette("bright",12)
color_order = [1,2,3,4,6,9,10]
colors = [palette_bright[i] for i in color_order]
#%%
# culmative coverage plot; methRate_df_filtered or methRate_df_filtered
plt.axvline(0.1, linestyle='--', color='black')
sns.ecdfplot(x='value', data=methRate_df, hue='group',#, weights='cov',
             palette=colors)
plt.savefig("Culm_rate_all.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()
plt.close()
#%%
plt.axvline(0.1, linestyle='--', color='black')
sns.ecdfplot(x='value', data=methRate_df_filtered, hue='group',#, weights='cov',
             palette=colors)
plt.xlim([0,1])
#plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))

plt.savefig("Culm_rate_all_filtered.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()
plt.close()

