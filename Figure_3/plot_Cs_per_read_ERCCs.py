import os, glob
import json
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

#%%
file = glob.glob("**/ERCC_C_Count_hist_nonDedup.json", recursive=True)[0]
with open(file, 'r') as data:
    C_counts = json.loads(data.read())


bin_list = [1, 2, 5, 10, 20, 150]

colnames = ['1', '2-4', '5-9', '10-19', '20+']
df = pd.DataFrame.from_dict(C_counts, orient='index')
df = df.sort_index()

bin_perc = df.div(df.sum(axis=1), axis=0)
bin_perc = bin_perc.sort_index()

#bin_perc = bin_perc.reset_index()

#%%
with open(file, 'r') as data:
    C_counts_all = json.loads(data.read())

df_all = pd.DataFrame.from_dict(C_counts_all, orient='index')
df_all = df_all.sort_index()

bin_perc = df_all.div(df_all.sum(axis=1), axis=0)
bin_perc['order'] = bin_perc.index
index_order = ['G5_1E_RNA_seq', 'G5_2E_RNA_seq', 'G1_1A', 'G1_2A', 'G2_1B', 'G2_2B',
               'G3_1C', 'G3_2C', 'G4_1D', 'G4_2D']
bin_perc['order'] = pd.Categorical(bin_perc.order, index_order)
bin_perc = bin_perc.sort_values('order')


#%%
fig, (ax) = plt.subplots(1,figsize=(4,4))

bin_perc.plot(kind='bar', stacked=True, color=['red','blue','green','orange','purple'],
              linewidth=1, edgecolor='black', ax=ax)
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
plt.savefig("C_content_bins_nonDedup.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()
plt.close()
