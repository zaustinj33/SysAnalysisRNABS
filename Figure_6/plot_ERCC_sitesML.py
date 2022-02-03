import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
#%%
#import sites info
total_df = pd.read_csv("Chapter_5\\ERCC_merged_in_SRR77all.csv", low_memory=False)  # file
#%%
# Plot methylation level
ERCC_df = total_df[total_df.filter(like='methRate').columns]

sns.stripplot(data=ERCC_df, order=['methRate_SRR8170377_Genome10xCall_1',
                                   'methRate_SRR8170377_Genome10xCall_2',
                                   'methRate_SRR8170377_Genome10xCall_3',
                                   'methRate_SRR8170377_Genome10xCall_5',
                                   'methRate_SRR8170377_Genome10xCall_8',
                                   'methRate_SRR8170377_Genome10xCall_15',
                                   'methRate_SRR8170377_Genome10xCall_all'])
plt.xticks(rotation=90)
#plt.savefig("ERCCsites_ML_SRR.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()

#%%

ERCC_df = total_df[total_df.filter(like='count').columns]

# Plot C-counts
sns.boxplot(data=ERCC_df, order=['C_count_SRR8170377_Genome10xCall_1',
                                 'C_count_SRR8170377_Genome10xCall_2',
                                 'C_count_SRR8170377_Genome10xCall_3',
                                 'C_count_SRR8170377_Genome10xCall_5',
                                 'C_count_SRR8170377_Genome10xCall_8',
                                 'C_count_SRR8170377_Genome10xCall_15',
                                 'C_count_SRR8170377_Genome10xCall_all'])

plt.xticks(rotation=90)
#plt.savefig("ERCCsites_counts_SRR.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()

#%%

# Plot number of unique sites
count_df = pd.DataFrame(ERCC_df.count()).T
#count_df = count_df
sns.barplot(data=count_df, order=['C_count_SRR8170377_Genome10xCall_1',
                                  'C_count_SRR8170377_Genome10xCall_2',
                                  'C_count_SRR8170377_Genome10xCall_3',
                                  'C_count_SRR8170377_Genome10xCall_5',
                                  'C_count_SRR8170377_Genome10xCall_8',
                                  'C_count_SRR8170377_Genome10xCall_15',
                                  'C_count_SRR8170377_Genome10xCall_all'])
plt.xticks(rotation=90)
plt.savefig("ERCCsites_siteCounts_SRR.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()
#%%

# Plot number of unique ERCCs per cutoff
count_df = total_df.copy()
test = count_df.groupby(['chrom']).mean()
test = test[test.filter(like='count').columns]
test_counts = pd.DataFrame(test.count()).T

sns.barplot(data=test_counts, order=['C_count_SRR8170377_Genome10xCall_1',
                                     'C_count_SRR8170377_Genome10xCall_2',
                                     'C_count_SRR8170377_Genome10xCall_3',
                                     'C_count_SRR8170377_Genome10xCall_5',
                                     'C_count_SRR8170377_Genome10xCall_8',
                                     'C_count_SRR8170377_Genome10xCall_15',
                                     'C_count_SRR8170377_Genome10xCall_all'])
plt.xticks(rotation=90)
plt.savefig("ERCCsites_ERCCs_perCutoff_SRR.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()

#%%
#count_df = count_df
sns.barplot(data=count_df, order=['C_count_SRR8170377_Genome10xCall_1',
                                  'C_count_SRR8170377_Genome10xCall_2',
                                  'C_count_SRR8170377_Genome10xCall_3',
                                  'C_count_SRR8170377_Genome10xCall_5',
                                  'C_count_SRR8170377_Genome10xCall_8',
                                  'C_count_SRR8170377_Genome10xCall_15',
                                  'C_count_SRR8170377_Genome10xCall_all'])
plt.xticks(rotation=90)
plt.savefig("ERCCsites_siteCounts_SRR.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()


#%%

# Plot number of unique ERCCs
ERCC_df = total_df.copy()
ERCC_df.index = ERCC_df['group']
ERCC_df = ERCC_df[ERCC_df.filter(like='all_all').columns]  # all sites

#ERCC_df = ERCC_df[ERCC_df.filter(like='1_Cutoff').columns].dropna()  # 1C cutoff
ERCC_df['group'] = ERCC_df.index
countERCC_df = pd.DataFrame(ERCC_df['group'].str.split().str[0].value_counts())
order_x = list(countERCC_df.index)  # define set for cutoff plots
y_max = max(pd.to_numeric(countERCC_df['group']))+10 # define ymax for cutoff plots

sns.barplot(data=countERCC_df.T, order=order_x)
plt.ylim(0,y_max)
plt.xticks(rotation=90)
plt.savefig("ERCCsites_ERCCcount_allC_SRR.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()
#%%

# Plot number of unique ERCCs by cutoff
ERCC_df = total_df.copy()
ERCC_df.index = ERCC_df['group']

ERCC_df = ERCC_df[ERCC_df.filter(like='_1').columns].dropna(0)  # 1C cutoff
ERCC_df['group'] = ERCC_df.index
countERCC_df = pd.DataFrame(ERCC_df['group'].str.split().str[0].value_counts()).T

# add back 0 counts
for ERCC in order_x:
    if ERCC not in countERCC_df.columns:
        countERCC_df[ERCC] = 0

sns.barplot(data=countERCC_df, order=order_x)
plt.xticks(rotation=90)
plt.ylim(0,y_max)
plt.savefig("ERCCsites_ERCCcount_1C_SRR.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()