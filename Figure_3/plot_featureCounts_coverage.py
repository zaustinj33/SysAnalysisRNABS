import glob, os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats

#%%

def get_FC_toDF(ext):
    """
    Calculates the FPKM or CPM of RNAseq counts files to compare deduplicated and non-deduplicated files.

    Args:
        ext: Extension (suffix) of featureCounts files
    Return:
         Dictionary of counts by geneID
    """
    colnames = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'cov']
    df_list = []
    name_list = []
    for file in glob.glob("**/*"+ext, recursive=True):
        print(file)
        name = os.path.basename(file).split(ext)[0]
        print(name)
        name_list.append(name)
        data_file = pd.read_csv(file, sep="\t", index_col=False, names=colnames, skiprows=2).iloc[:,[0,2,5,6]]
        #data_file['normCov'] = data_file['cov'].div((data_file['cov'].sum()*data_file['Length'])) * 1e6 # FPKM
        data_file['normCov'] = data_file['cov'].div((data_file['cov'].sum()) * 1e6) # CPM
        fpkm_file = data_file.iloc[:,[0,3,4]]

        fpkm_file = fpkm_file[fpkm_file.iloc[:,2] > 0]  # FPKM > 0
        df_list.append(fpkm_file)

    #result = pd.concat(df_list, keys=name_list)
    return dict(zip(name_list, df_list))
nonDedup = get_FC_toDF("_genomeCounts.txt")
Dedup = get_FC_toDF("_genomeCountsDedup.txt")


#%%
fpkm_list = []
mergedName_list = []
for k in Dedup:
    fpkm_list.append(Dedup[k].merge(nonDedup[k], on='Geneid', suffixes=['Dedup', 'nonDedup']))
    mergedName_list.append(k)

result = pd.concat(dict(zip(mergedName_list, fpkm_list)), axis=0).reset_index()
result['group'] = result['level_0'].str.replace(r'_.*$', '')  # add group level
#%%
result['rep'] = result['level_0'].str.replace('^.*_','')  # add group level
result['rep'] = result['rep'].str.replace(r'[A-Z]+','')
#%%

result['change'] = abs(result['normCovnonDedup'] - result['normCovDedup']) / result['normCovnonDedup']
result['change_bin'] = np.where(result['change'] > 1, 'High Change', 'Low Change') # LFC > 1

result_ERCC = result[result['Geneid'].str.contains("ERCC")]
#%%
names_colors = sns.color_palette("bright").as_hex()
names_colors.pop(0)

b1 = sns.boxplot(data=result, x='group', y='covDedup', palette=names_colors, orient='v')
b1.set(yscale='log')
#plt.savefig("Boxplot_CovDedup.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()

#%%
b2 = sns.boxplot(data=result_ERCC, x='Geneid', y='covnonDedup', palette=names_colors, orient='v')
b2.set(yscale='log')
#plt.savefig("Boxplot_CovNonDedup_ERCC.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()
#%%
g = sns.FacetGrid(result, col='group', hue='change_bin', palette=['black', 'gold'])
g.map(sns.scatterplot, 'normCovnonDedup', 'normCovDedup', size=1, alpha=0.5)#, hue='group', palette='bright',alpha=0.1)

# draw y=x
"""
x0, x1 = p.get_xlim()
y0, y1 = p.get_ylim()
lims = [max(x0, y0), min(x1, y1)]
p.plot(lims, lims, '-r')
"""
g.set(xscale='log')
g.set(yscale='log')

#plt.savefig("featureCounts_coverage.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()


