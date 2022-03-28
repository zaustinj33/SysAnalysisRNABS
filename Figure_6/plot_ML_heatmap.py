import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd
from scipy.spatial import distance
from scipy.cluster import hierarchy
from sklearn.cluster import AgglomerativeClustering

os.chdir('Chapter_5')

# %%
# import all sites
conditions = ["G", "_MF", "pMF"]  # ,"SRR"]

total_file = "allm5C_libraries_filteredDepthAnno.csv"
total_df = pd.read_csv(total_file, low_memory=False)  # file
total_df = total_df.sort_values(by=['position'])  # sort
# %%
names = ['G1', 'G2', 'G3', 'G4', 'MF_rep1', 'MF_rep2', 'pMF_rep1', 'pMF_rep2', 'rep1',
         'rep2', 'rep3', 'rep4']
# Aggregate methylation level for each condition
total_df.index = total_df['group']
cov_df = total_df.filter(regex='cov')
count_df = total_df.filter(regex='count')

cov_dict = {}
count_dict = {}
for name in conditions:
    cov_dict[name] = cov_df.filter(regex=name).sum(axis=1)
    count_dict[name] = count_df.filter(regex=name).sum(axis=1)

ML_dict = {}

for i, j in cov_dict.items():
    ML_dict[i] = count_dict[i].divide(j, fill_value=0)

result_df = pd.DataFrame(ML_dict)

# result_df.dropna(axis=0, inplace=True, subset=['SRR','_MF','pMF'])
# result_df.replace(np.nan, 0, inplace=True)
# result_df.replace(0, np.nan, inplace=True)

result_df = result_df[(result_df['G'] > 0.1) | (result_df['_MF'] > 0.1) |
                      (result_df['pMF'] > 0.1)]  # | (result_df['SRR'] > 0.1)]

result_df.dropna(axis=0, inplace=True)
test = total_df[total_df['group'].isin(result_df.index)]
# test.to_csv("AllConditionOverlap_methylationLevel.csv")

# %%
result_df_ML = total_df.filter(regex="methRate")
result_df_ML.replace(np.nan, 0, inplace=True)

cov_df.columns = names
count_df.columns = names
# %%
from matplotlib.colors import LinearSegmentedColormap

boundaries = [0.0, 0.05, 0.1, 0.2, 0.4, 0.6, 1.0]
hex_colors = sns.color_palette("RdYlBu_r", n_colors=len(boundaries) * 2).as_hex()
hex_colors = [hex_colors[i] for i in range(0, len(hex_colors), 2)]

colors = list(zip(boundaries, hex_colors))

custom_color_map = LinearSegmentedColormap.from_list(
    name="cus",
    colors=colors,
)
# %%
# Define clusters
correlations_array = np.asarray(result_df)

row_linkage = hierarchy.linkage(
    distance.pdist(correlations_array), method='ward')

col_linkage = hierarchy.linkage(
    distance.pdist(correlations_array.T), method='ward')

model = AgglomerativeClustering(n_clusters=8, affinity='euclidean', linkage='ward')
model = model.fit_predict(correlations_array)

# %%
lut = dict(zip(set(model), ['red', 'blue', 'green', 'orange', 'purple', 'pink', 'black', 'grey']))
row_colors = pd.DataFrame(model)[0].map(lut)

cg = sns.clustermap(result_df.reset_index(drop=True), row_linkage=row_linkage, col_linkage=col_linkage,
                    cmap=custom_color_map,
                    row_colors=row_colors, figsize=(5, 5), yticklabels=False, col_cluster=False,
                    robust=True, method='ward')  # , row_cluster=False)  # z_score=0,
cg.ax_row_dendrogram.set_visible(False)
# plt.savefig("ML_conditions_clusteringHeatmapDepth.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()

plt.close()

# %%
merge_df = result_df
merge_df['cluster'] = model
merge_df['group'] = result_df.index
merge_df.reset_index(drop=True)
cluster_df = pd.merge(merge_df.rename_axis(None), total_df.rename_axis(None), on='group')

cluster_gene_list = (cluster_df['gene_name'][cluster_df['cluster'] == 5]).unique()

cluster_file = open("Total_cluster_genes.txt", "w")
for i in cluster_gene_list:
    cluster_file.write(i + '\n')
cluster_file.close()

# %%
from scipy.stats import zscore

# write correlation matrix (z-score)

zscore_vals = result_df.apply(zscore, axis=1)

# %%
from scipy import stats


# BH t-test
def BH_test(set1, set2):
    # subset tests by relevant sites identified by 04a_OverlapDotplot.R
    master_set = pd.read_csv('Dotplot_' + set1 + set2 + '_table.csv')
    master_set = master_set.dropna(subset=['ML_1', 'ML_2']).reset_index()
    count_set = {set1: master_set['C_count_' + set1], set2: master_set['C_count_' + set2]}
    cov_set = {set1: master_set['cov_' + set1], set2: master_set['cov_' + set2]}

    pvals = []
    p_adj = []
    try:
        len(count_set[set1]) == len(cov_set[set1])
    except:
        print('data is not same size')
    for i in range(len(count_set[set1])):
        cont_table = pd.DataFrame({set1: [count_set[set1][i], cov_set[set1][i]],
                                   set2: [count_set[set2][i], cov_set[set2][i]]})

        odds, pvalue = stats.fisher_exact(cont_table)
        pvals.append(pvalue)
    pvals_sorted = sorted(pvals, key=float)  # sorted pvalues
    master_set['pval'] = pvals
    master_set = master_set.sort_values('pval', ascending=True)

    rank = 1
    for p in pvals_sorted:
        fdr_pval = p * len(pvals_sorted) / rank
        rank += 1
        p_adj.append(fdr_pval)

    master_set['BH'] = p_adj
    master_set['shape'] = np.where(master_set['BH'] <= 0.01, 'sig', 'non-sig')

    return master_set


test_BH = pd.DataFrame(BH_test('G3', 'G4'))
# %%
rcParams['figure.figsize'] = 3, 3
markers = {"sig": "X", "non-sig": "o"}
palette = ['blue']
# ax = sns.scatterplot(data=test_BH[test_BH['BH'] > 0.01], x='ML_1', y='ML_2', style = 'shape',
#                    markers=markers, s=25)

sns.scatterplot(data=test_BH, x='ML_1', y='ML_2', style='shape', hue='shape', palette = palette,
                markers=markers, s=25)
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.legend([], frameon=False)
plt.savefig("G3G4_DMS.png",bbox_inches='tight', dpi=400, transparent=True)
plt.show()

# %%

# Correlation matix of samples
from scipy.spatial import distance
from scipy.cluster import hierarchy

correlations = result_df.corr()
correlations_array = np.asarray(result_df.corr())

row_linkage = hierarchy.linkage(
    distance.pdist(correlations_array), method='average')

col_linkage = hierarchy.linkage(
    distance.pdist(correlations_array.T), method='average')

sns.clustermap(correlations, row_linkage=col_linkage, col_linkage=row_linkage, method="average",
               figsize=(5, 10))
plt.show()

# %%
from matplotlib.colors import LinearSegmentedColormap

boundaries = [0.0, 0.05, 0.1, 0.2, 0.4, 0.6, 1.0]
hex_colors = sns.color_palette("RdBu_r", n_colors=len(boundaries) * 2 + 2).as_hex()
hex_colors = [hex_colors[i] for i in range(0, len(hex_colors), 2)]

colors = list(zip(boundaries, hex_colors))

custom_color_map = LinearSegmentedColormap.from_list(
    name="cus",
    colors=colors,
)

cg = sns.clustermap(result_df, annot=False, cmap=custom_color_map, dendrogram_ratio=(.1, .2),
                    figsize=(5, 5), yticklabels=False)  # z_score=0,
cg.ax_row_dendrogram.set_visible(False)
plt.savefig("ML_conditions_clusteringHeatmapCcutoffDepth_noSRR.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()

plt.close()
