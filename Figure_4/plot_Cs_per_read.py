import os
import json
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


#%%
"""
Plot the output of Preprocessing/Read Processing/record_Cs_perRead.py
Python file can be modified to record different bins, files, etc
"""

## Change input to file of interest (ERCC_Cs, noDedup_ERCCs, etc)
with open((os.getcwd()+'\\C_Count_hist_nonMT.json'), 'r') as data:
    C_counts = json.loads(data.read())


bin_list = [1, 2, 5, 10, 20, 150]

colnames = ['1', '2-4', '5-9', '10-19', '20+']
df = pd.DataFrame.from_dict(C_counts, orient='index')
df = df.sort_index()
#%%
bin_perc = df.div(df.sum(axis=1), axis=0)
bin_perc.index = bin_perc.index.str.replace(r'_meRanGh.*$', '')
bin_perc['order'] = bin_perc.index

index_order = ['G5_1E_RNA_seq', 'G5_2E_RNA_seq', 'G1_1A', 'G1_2A', 'G2_1B', 'G2_2B',
               'G3_1C', 'G3_2C', 'G4_1D', 'G4_2D', 'MFRNAseq_rep1', 'MFRNAseq_rep2', 'MF_rep1',
               'MF_rep2', 'pMFRNAseq_rep1', 'pMFRNAseq_rep2', 'pMF_rep1', 'pMF_rep2',
               'SRR8170377', 'SRR8170378', 'SRR8170379', 'SRR8170380']
bin_perc['order'] = pd.Categorical(bin_perc.order, index_order)

bin_perc = bin_perc.sort_values('order')
#bin_perc.sort_values(index_order)
#%%
fig, (ax) = plt.subplots(1,figsize=(4,4))

# Drop 20+ for zoomed in view
#bin_perc_small = bin_perc.drop(columns=['20+'])

bin_perc.plot(kind='bar', stacked=True, color=['grey','red','blue','green','orange','purple'],
              edgecolor='black', width=0.5, linewidth=1, ax=ax)
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))


plt.savefig("C_content_bins_ERCC_Dedup.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()
plt.close()

