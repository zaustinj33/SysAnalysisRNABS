## plot read annotation distrubition from RSeQC file

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import os, glob
import numpy as np

#%%
df_list = []
name_list = []
for file in glob.glob("**/*_DistributionSummary.txt", recursive=True):  # remove '**/' if files are local to cwd
    name = os.path.basename(file).split("_DistributionSummary.txt")[0]
    print(name)
    df = pd.read_csv(file, sep='\\s+', skiprows=4)[0:10]  # import weird tsv
    df.index=df['Group']
    df = df.iloc[: , 1:]  # skip weird rows
    df['perc_tag'] = df['Tag_count']/df['Tag_count'].sum()

    df_intergenic = df[df.index.str.contains("kb")].sum()  # add intergenic as sum of upstream of TSS and down of TES
    df_intergenic.rename("intergenic")

    df = df.append(df_intergenic, ignore_index=True)
    df.index = ['Exon', '5\'UTR', '3\'UTR', 'Intron', 'TSS_up_1kb','TSS_up_5kb','TSS_up_10kb',
                'TES_down_1kb', 'TES_down_5kb', 'TES_down_10kb', 'Intergenic']

    df_list.append(df)
    name_list.append(name)

result = pd.concat(df_list, keys=name_list)

#%%

result.index.names = ['name', 'anno']
# melt dataframe for better plotting
melt_df = result.stack(level=0).reset_index(level=2, drop=False).reset_index()
melt_df['group'] = melt_df['name'].str.replace(r'817.*$|_.*$', '')  # add levels
melt_df = melt_df[melt_df['level_2'] == 'perc_tag']
plot_df = melt_df[~melt_df['anno'].str.contains('kb')]

order_x = ['G5', 'G1', 'G2', 'G3', 'G4', 'MFRNAseq', 'MF', 'pMFRNAseq', 'pMF', 'SRR']

#%%
"""
test_values = df.iloc[0:,2]
test_labels = ['Exon', '5\'UTR', '3\'UTR', 'Intron', 'Intergenic']

# Pie chart for single sample
#define Seaborn color palette to use
colors = sns.color_palette('bright')[0:5]

#create pie chart
plt.pie(test_values, labels=test_labels,
        colors=colors, autopct='%.0f%%')
plt.show()
"""
#%%
# Barchart for all samples
sns.despine()
sns.set_style("ticks")

colors = sns.color_palette('bright')[0:5]
p = sns.catplot(x='group', y=0, hue='anno', data=plot_df, size=4, aspect=2,
                edgecolor="black", linewidth=2, palette=colors, kind='bar', order=order_x)

p.set(yticks=np.arange(0, 1.1, 0.5))
plt.xticks(rotation=90)
plt.legend([],[], frameon=False)
plt.xlabel(None)

plt.savefig('Perc_IntronExon_fromBAM.png', bbox_inches='tight', dpi=400, transparent=True)
plt.show()
plt.close()

#%%

plot_df[0][(plot_df['anno'].str.contains('Exon')) & plot_df['group'].str.contains('RNAseq|G5')].mean() + \
plot_df[0][(plot_df['anno'].str.contains('5\'UTR')) & plot_df['group'].str.contains('RNAseq|G5')].mean() + \
plot_df[0][(plot_df['anno'].str.contains('3\'UTR')) & plot_df['group'].str.contains('RNAseq|G5')].mean()

