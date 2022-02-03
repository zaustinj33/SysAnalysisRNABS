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
conditions = ["G","MF","pMF","SRR"]

total_file = "Chapter_5\\allm5C_libraries_filteredDepthAnno.csv"
total_df = pd.read_csv(total_file, low_memory=False)  # file

#%%

# all sites' distribution
binLocDist = total_df['binLoc'].value_counts().reset_index()

# Create data interpolation
x_new = np.linspace(1, 24, 100)
a_BSpline = interpolate.make_interp_spline(binLocDist['index'].sort_values(), binLocDist['binLoc'])
y_new = a_BSpline(x_new)

interp_df = pd.DataFrame(np.row_stack([x_new, y_new])).T
#%%

# distribution for overlapped sites
overlapped_site_bins = {}
for set in conditions:
    print(set)
    for file in glob.glob("**/"+set+"_allOverlapDepthAnno.csv", recursive=True):
        df_set = pd.read_csv(file, low_memory=False)
        overlapped_site_bins[set] = df_set['binLoc'].value_counts()

#%%

overlapped_site_df = pd.DataFrame(overlapped_site_bins)
overlapped_site_df.replace(np.nan, 0, inplace=True)
overlapped_site_df = overlapped_site_df / overlapped_site_df.sum()
overlapped_site_df = overlapped_site_df.reset_index()

#%%
# Create data interpolation
x_new = np.linspace(1, 24, 50)
interp_list = []
for column in overlapped_site_df.columns:
    print(column)
    a_BSpline = interpolate.make_interp_spline(overlapped_site_df['index'].sort_values(), overlapped_site_df[column])
    y_new = a_BSpline(x_new)
    interp_list.append(y_new)
#%%
interp_df = pd.DataFrame(interp_list).T
interp_df.columns = overlapped_site_df.columns
interp_df = interp_df.drop(['index'], axis=1)

#%%
site_bin_melt = overlapped_site_df.stack().reset_index()
#%%
site_bin_melt['name'] = np.where(site_bin_melt['level_1'] == 'G', 'MT',
                                 np.where(site_bin_melt['level_1'] == 'MF', 'Total RNA',
                                 np.where(site_bin_melt['level_1'] == 'pMF', 'Polysome RNA',
                                 np.where(site_bin_melt['level_1'] == 'SRR', 'Huang', 'Index'))))
#%%
rcParams['figure.figsize'] = 5,3

#sns.scatterplot(data=binLocDist.reset_index(), x='index', y='binLoc', color='blue')
plt.plot(interp_df[0], interp_df[1], color='red', linewidth=3)

#sns.lineplot(data=site_bin_melt, x='level_1', y=0, hue='name')
plt.axvline(5, 1,0, linestyle='--',linewidth=2, color='black')
plt.axvline(17, 1,0, linestyle='--',linewidth=2, color='black')
plt.xlim(0,25)
plt.ylim(0,max(binLocDist['binLoc']+100))
plt.legend([], frameon=False)
plt.savefig("site_distribution.png",bbox_inches='tight', dpi=400, transparent=True)
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
conditions = ["G","_MF","pMF","SRR"]

total_file = "Chapter_5\\allm5C_libraries_filteredDepthAnno.csv"
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
