import glob
import pandas as pd
import os, re
import numpy as np
import matplotlib.pyplot as plt

#%%
"""
Plots binned 1) coverage 2) C-counts for all unfiltered meRanCall files
"""
headers = ['#SeqID','refPos', 'refStrand', 'refBase', 'cov', 'C_count', 'methRate', 'mut_count',
           'mutRate', 'CalledBase', 'CB_count',	'state', '95_CI_lower', '95_CI_upper', 'p-value_mState',
           'p-value_mRate', 'score', 'seqContext', 'geneName', 'candidateName', 'gene']
count_dict = {}
coverage_dict = {}
ML_dict = {}

for file in glob.glob("**/*_Genome10xCall_annotate.txt", recursive=True):  #*_Genome10xCall*_annotate.txt", recursive=True):
    name = re.sub("_Genome10xCall_annotate.txt", "", os.path.basename(file))
    print(file)
    df = pd.read_csv(file, sep='\t', low_memory=False)  # file
    df.columns = headers
    df = df[df['state'] == 'M']

    #1) basic filter
    df['binned_cov'] = np.where(df['cov'] < 20, '10-19',
                                np.where(df['cov'] < 100, '20-99',
                                np.where(df['cov'] < 1000, '100-999',
                                np.where(df['cov'] < 10000, '1,000-9,999', '10,000+'))))

    df['binned_count'] = np.where(df['C_count'] < 3, '1-2',
                            np.where(df['C_count'] < 10, '3-9',
                            np.where(df['C_count'] < 20, '10-19',
                            np.where(df['C_count'] < 100, '20-99','100+'))))

    df['binned_ML'] = np.where(df['methRate'] < 0.1, '< 0.1',
                        np.where(df['methRate'] < 0.2, '< 0.2',
                        np.where(df['methRate'] < 0.4, '< 0.4',
                        np.where(df['methRate'] < 1, '< 1','Error'))))

    coverage_dict[name] = df['binned_cov'].value_counts()
    count_dict[name] = df['binned_count'].value_counts()
    ML_dict[name] = df['binned_ML'].value_counts()

#%%
cov_df = pd.DataFrame(coverage_dict).T
count_df = pd.DataFrame(count_dict).T
ML_df = pd.DataFrame(ML_dict).T

# Sort index
bin_perc = cov_df.div(cov_df.sum(axis=1), axis=0)
bin_perc.index = bin_perc.index.str.replace(r'_meRanGh.*$', '')
bin_perc['order'] = bin_perc.index

#index_order = ['G5_1E_RNA_seq', 'G5_2E_RNA_seq', 'G1_1A', 'G1_2A', 'G2_1B', 'G2_2B',
 #              'G3_1C', 'G3_2C', 'G4_1D', 'G4_2D', 'MFRNAseq_rep1', 'MFRNAseq_rep2', 'MF_rep1',
  #             'MF_rep2', 'pMFRNAseq_rep1', 'pMFRNAseq_rep2', 'pMF_rep1', 'pMF_rep2',
   #            'SRR8170377', 'SRR8170378', 'SRR8170379', 'SRR8170380']

index_order = ['G1', 'G2', 'G3', 'G4']
bin_perc['order'] = pd.Categorical(bin_perc.order, index_order)

bin_perc = bin_perc.sort_values('order')

# Sort legend
cov_order = ['10-19','20-99','100-999','1,000-9,999','10,000+']
bin_perc.columns = pd.CategoricalIndex(bin_perc.columns.values,
                                       ordered=True,
                                       categories=cov_order)
bin_perc = bin_perc.sort_index(axis=1)
#%%
cov_bin_colors = ['#c20000', '#ff3d00', '#ff7600', '#ffb74d', '#ffd900']
plt.rcParams["figure.figsize"] = (20, 5)

bin_perc.plot(kind='bar', stacked=True, color=cov_bin_colors,
              edgecolor='black', width=0.5, linewidth=2)
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
plt.savefig("Coverage_bins_Ccutoff.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()
plt.close()

#%%
count_bin_colors = ['#d50000', '#ff5131', '#aa00ff', '#2962ff', '#64dd17']
bin_perc = count_df.div(count_df.sum(axis=1), axis=0)
bin_perc.index = bin_perc.index.str.replace(r'_meRanGh.*$', '')

# Sort index
bin_perc['order'] = bin_perc.index
bin_perc['order'] = pd.Categorical(bin_perc.order, index_order)
bin_perc = bin_perc.sort_values('order')

# Sort legend
count_order = ['1-2','3-9','10-19','20-99','100+']
bin_perc.columns = pd.CategoricalIndex(bin_perc.columns.values,
                                 ordered=True,
                                 categories=count_order)
bin_perc = bin_perc.sort_index(axis=1)

bin_perc.plot(kind='bar', stacked=True, color=count_bin_colors,
              edgecolor='black', width=0.5, linewidth=2)
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
plt.savefig("Count_bins_Cutoff.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()
plt.close()

#%%
ML_bin_colors = ['#d50000', '#ff5131', '#aa00ff', '#2962ff', '#64dd17']
bin_perc = ML_df.div(count_df.sum(axis=1), axis=0)

bin_perc.index = bin_perc.index.str.replace(r'_meRanGh.*$', '')

# Sort index
bin_perc['order'] = bin_perc.index
bin_perc['order'] = pd.Categorical(bin_perc.order, index_order)
bin_perc = bin_perc.sort_values('order')

# Sort legend
count_order = ['< 0.1', '< 0.2', '< 0.4', '< 1']
bin_perc.columns = pd.CategoricalIndex(bin_perc.columns.values,
                                       ordered=True,
                                       categories=count_order)
bin_perc = bin_perc.sort_index(axis=1)

bin_perc.plot(kind='bar', stacked=True, color=count_bin_colors,
              edgecolor='black', width=0.5, linewidth=2)
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
plt.savefig("ML_bins.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()
plt.close()
