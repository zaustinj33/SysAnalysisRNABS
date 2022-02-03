import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import glob, os
#from calc_UMI_errorRate import *
import numpy as np
import pickle

#%%
#from Chapter_3.calc_UMI_errorRate import subset_bam  ## not working

counts_list = {}
coverage_list = {}
name_list = []
for file in glob.glob("**/*_meRanGh_genomeMapTrim_dedupGroup.bam_ERCC.sam", recursive=True):
    print(file)
    name = os.path.basename(file).split("_meRanGh_genomeMapTrim_dedupGroup.bam_ERCC.sam")[0]
    name_list.append(name)
    counts_file, coverage_file = countCs_from_SAM(file)
    #subset_bam(file)

    counts_list[name] = pd.Series(counts_file.values())
    coverage_list[name] = pd.Series(coverage_file.values())
    name_list.append(name)
#%%
coverage_df = pd.DataFrame(coverage_list)
coverage_melt = coverage_df.melt()
bins = [0,1,4,9,19,1000]
labels = ['1','2-4','5-9','10-19','20+']
coverage_melt['binned'] = pd.cut(coverage_melt['value'], bins, labels=labels)

#%%
coverage_binned = coverage_melt.pivot(values='binned', columns='variable')
coverage_binned = coverage_binned.apply(pd.value_counts, axis=0).T
coverage_binned.index.name = 'sample'
coverage_binned.columns=coverage_binned.columns.astype('str')
coverage_binned = coverage_binned.reset_index()
coverage_binned_melt = pd.melt(coverage_binned, id_vars='sample')
coverage_binned_melt['group'] = coverage_binned_melt['sample'].str.replace(r'_.*$', '')
#%%
#sns.stripplot(data=coverage_df[coverage_df > 0], palette='Paired')
sns.boxplot(data=coverage_df[coverage_df > 1], palette='Paired')
plt.show()

#%%
p = sns.catplot(data=coverage_binned_melt, x='group', y='value', hue='variable', size=4, aspect=2,
                edgecolor="black", linewidth=2, palette='bright', kind='bar')
plt.xticks(rotation=90)
plt.yscale('log')
sns.despine()
sns.set_style("ticks")

#plt.savefig("UMI_group_coverage.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()
#%%
# m5C rate dict
rate_dict = []
counts_dict = []
for k, v in counts_list.items():  # k = sample name
    #print(k)
    for i, j in v.items():  # i = records on read
        #print(j)
        for x,y in j.items():  # x = position on read, y = nucleotide coverage at position in group
            #print(y)
            if ((coverage_list[k][i] >= 1) & (y >=1)):  ## IMPORTANT GROUP COVERAGE, nuc count, CUTOFF ##
                # (name of sample, coverage at position, count, m5C rate)
                rate_dict.append(tuple([k, coverage_list[k][i], int(y), (int(y) / coverage_list[k][i])]))
#%%
rate_df = pd.DataFrame.from_records(rate_dict, columns=['name', 'coverage', 'count', 'rate'])
rate_df['group'] = rate_df['name'].str.replace(r'817.*$|_.*$', '')  # add levels

rate_df['conv_group'] = np.where((rate_df['rate'] == 1) | (rate_df['rate'] == 0), 'Concordant reads',
                                 'Discordant reads')
rate_df = rate_df[rate_df['coverage'] > 1]
#%%
# save dataframe as table
cwd = os.getcwd()
path = cwd + "/T_discordance_stats_plus2Cov.pkl"
rate_df.to_pickle(path)
#%%
## "C_discordance_stats_allCov.pkl" ## "C_discordance_stats.pkl"
rate_df = pickle.load(open("T_discordance_stats_plus2Cov.pkl", 'rb'))

#%%
conv_rate_bin_df = rate_df.groupby(['name','group'])['conv_group'].value_counts().reset_index(name='count')

#%%
order_x = ['G5', 'G1', 'G2', 'G3', 'G4']
# Conversion efficiency bargraph
sns.barplot(data=conv_rate_bin_df, x='group', y='count',
            hue='conv_group', order=order_x, linewidth=2, edgecolor='black')
plt.yscale('log')
plt.ylim(1, 10e7)
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
plt.xticks(rotation=90)

plt.savefig("Concord_Discord_plus2Cov_C.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()

#%%

# C_count histogram
sns.histplot(rate_df['count'][rate_df['name'] == 'G4_1D'])
plt.yscale('log')
plt.show()
#%%
# coverage/error rate curve
sns.scatterplot(data=rate_df[rate_df['name'] == 'G4_1D'], x='coverage', y='rate', hue='name', palette='bright')
plt.show()

#%%

# coverage of ERCC non-conversion sites
sns.boxplot(y=rate_df['coverage'], x=rate_df['name'])
plt.show()

#%%

# rates of ERCC non-conversion sites
sns.boxplot(y='rate', x='name', data=rate_df[(rate_df['rate'] != 1) & (rate_df['rate'] != 0)])
plt.show()


