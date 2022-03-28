import glob
import json, os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import interpolate
from matplotlib import rcParams

#%%

#import all sites
conditions = ["G","SRR"] #"MF","pMF",

total_file = "allm5C_libraries_filteredDepthAnno.csv"
total_df = pd.read_csv(total_file, low_memory=False)  # file

# %%
# all sites' distribution
binLocDist = total_df['binLoc'].value_counts().reset_index()

# %%

# distribution for overlapped sites
overlapped_site_bins = {}
for set in conditions:
    print(set)
    for file in glob.glob("**/" + set + "_allOverlapRepDepth.csv", recursive=True):
        df_set = pd.read_csv(file, low_memory=False)
        df_set_anno = total_df[total_df['group'].isin(df_set['group'])]
        overlapped_site_bins[set] = df_set_anno['binLoc'].value_counts()

# %%
overlapped_site_df = pd.DataFrame(overlapped_site_bins).T
overlapped_site_df.replace(np.nan, 0, inplace=True)
site_bin_melt = overlapped_site_df.stack().reset_index()
site_bin_melt = site_bin_melt[(site_bin_melt['level_1'] < 24) & (site_bin_melt['level_1'] >= 0)]

site_bin_interp = site_bin_melt.copy()
site_bin_interp['name'] = pd.Series(site_bin_interp['level_0'].str.match(r'^SRR.*'))
site_bin_interp['fraction'] = np.where(site_bin_interp['name'] == True, "SRR", "MT")
#%%
from scipy.interpolate import interp1d
site_bin_interp_test = site_bin_interp[site_bin_interp['fraction'] == 'MT']
# Interpolate points
test = site_bin_interp_test.groupby('level_1').mean()
f_cubic = interp1d(test.index, test[0], kind='cubic')
xnew = np.linspace(1, 23, num=300, endpoint=True)
#%%
plt.plot(xnew, f_cubic(xnew), '-', label='cubic', color='blue')
plt.axvline(5, 1, 0, linestyle='--', linewidth=2, color='black')
plt.axvline(17, 1, 0, linestyle='--', linewidth=2, color='black')
plt.xlim(0, 24)
plt.ylim(0,12)
plt.yticks(np.arange(0, 12+1, 12/2), labels=[])
plt.xticks([])
plt.savefig("site_distribution_condition_MT.png",bbox_inches='tight', dpi=200, transparent=True)
plt.show()

# %%
rcParams['figure.figsize'] = 5, 3

#sns.scatterplot(data=binLocDist.reset_index(), x='index', y='binLoc', color='blue')
plt.plot(interp_df[0], interp_df[1], color='red', linewidth=3)

# sns.lineplot(data=site_bin_melt, x='level_1', y=0, hue='name')
plt.axvline(5, 1, 0, linestyle='--', linewidth=2, color='black')
plt.axvline(17, 1, 0, linestyle='--', linewidth=2, color='black')
plt.xlim(0, 25)
#plt.ylim(0, max(binLocDist['binLoc'] + 100))
plt.legend([], frameon=False)
# plt.savefig("site_distribution_trueOverlap.png",bbox_inches='tight', dpi=400, transparent=True)
plt.show()

plt.close()

#%%

sns.lineplot(data=site_bin_melt, x='level_0', y=0, hue='name')
#plt.axvline(5, 1,0, linestyle='--')
#plt.axvline(17, 1,0, linestyle='--')
plt.xlim(0,21)
plt.ylim(0,0.2)
#plt.legend([], frameon=False)
plt.savefig("site_distribution_conditions.png",bbox_inches='tight', dpi=400, transparent=True)
plt.show()

plt.close()

#%%
#import all sites
conditions = ["G","SRR"] #"_MF","pMF",

total_file = "allm5C_libraries_filteredDepthAnno.csv"
total_df = pd.read_csv(total_file, low_memory=False)  # file
#%%

def stack_ML(df):
    ML_df = df.filter(regex='methRate').stack()
    return ML_df

UTR5_ML = stack_ML(total_df[total_df['type'] == '5UTR'])
UTR3_ML = stack_ML(total_df[total_df['type'] == '3UTR'])
CDS_ML = stack_ML(total_df[total_df['type'] == 'exon'])
#%%

#names=['5UTR','CDS','3UTR']
plot_dist_ML = pd.DataFrame([UTR5_ML,CDS_ML,UTR3_ML]).T

sns.violinplot(data=plot_dist_ML, linewidth=2, palette='husl',cut=True)
plt.ylim(-0.1,1.1)
plt.savefig("ML_siteDist_violinplot.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()
