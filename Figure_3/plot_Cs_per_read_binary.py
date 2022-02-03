import os
import json
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


#%%

## Change input to file of interest (ERCC_Cs, noDedup_ERCCs, etc)
with open((os.getcwd()+'\\C_Count_hist_binary_noDedup.json'), 'r') as data:
    C_counts = json.loads(data.read())

colnames = ['no C', 'with C']
df = pd.DataFrame.from_dict(C_counts, orient='index')
df = df.sort_index()

bin_perc = df.div(df.sum(axis=1), axis=0)
bin_perc.index = bin_perc.index.str.replace(r'_meRanGh.*$', '')
bin_perc['order'] = bin_perc.index

index_order = ['G1_1A', 'G1_2A', 'G2_1B', 'G2_2B',
               'G3_1C', 'G3_2C', 'G4_1D', 'G4_2D', 'MF_rep1',
               'MF_rep2', 'pMF_rep1', 'pMF_rep2',
               'SRR8170377', 'SRR8170378', 'SRR8170379', 'SRR8170380']
bin_perc['order'] = pd.Categorical(bin_perc.order, index_order)
bin_perc['group'] = bin_perc['order'].str.replace(r'817.*$|_.*$', '')

bin_perc = bin_perc.sort_values('order')
#bin_perc.sort_values(index_order)
#%%
pal = sns.color_palette('bright')
slices = [1,2,3,4,6,8,9]
colors = [pal.as_hex()[i] for i in slices]

fig, (ax) = plt.subplots(1,figsize=(3,4))

p = sns.barplot(data=bin_perc, x='group', y='no C', palette=colors, ax=ax,
                linewidth=1, edgecolor='black')
plt.xticks(rotation=90)
plt.ylim(0.9,1)
p.set(yticks=np.arange(0.8, 1.01, 0.1))
sns.despine()

plt.savefig("C_content_binary_noDedup.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()
plt.close()

