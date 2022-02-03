import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#%%
# Create ML of sites overlapped ([# replicates in set] + 1)
def calcML(df):
    ML_df = df[df.filter(like="methRate").columns].stack()
    return ML_df

#samples = ['G1' 'G2', 'G3', 'G4','MF_rep1','MF_rep2','SRR8170377', 'SRR8170378', 'SRR8170379','SRR8170380']
conditions = ["G","MF","pMF","SRR"]
#%%
#import all sites
total_file = "allm5C_libraries_filteredDepthAnno.csv"
total_df = pd.read_csv(total_file, low_memory=False)  # file

overlapped_site_ML = {}
exclusive_site_ML = {}
for set in conditions:
    print(set)
    """if set == "G":
        ext = "_allUnionDepth.csv"
    else:"""

    ext = "_allOverlapDepth.csv"
    #ext = "_allUnion.csv"
    for file in glob.glob("**/"+set+ext, recursive=True):
        print(file)
        df_set = pd.read_csv(file, low_memory=False)
        name = "_"+set
        df = df_set[df_set.filter(like=name).columns]
        df['group'] = df_set['group']
        df = df[df['group'].isin(total_df['group'])]

        # subset total site columns df by name
        df_exclusive = total_df[total_df.filter(like=name).columns]
        df_exclusive['group'] = total_df['group']
        # Remove sites
        df_exclusive = df_exclusive[~df_exclusive['group'].isin(df_set['group'])]

        print("subset overlapped")
        overlapped_site_ML[set] = calcML(df)
        print("subset exlusive")
        exclusive_site_ML[set] = calcML(df_exclusive)

overlapped_df = pd.DataFrame(overlapped_site_ML)
exclusive_df = pd.DataFrame(exclusive_site_ML)
#%%
flierprops = dict(marker='o', markerfacecolor='black', markersize=3,
                  linestyle='none', markeredgecolor='black')

sns.boxplot(data=overlapped_df, linewidth=2,flierprops=flierprops)
#sns.scatterplot(data=overlapped_df, color='black')
plt.ylim(-0.05,1.1)
plt.savefig("overlapped_ML_boxplot.png",bbox_inches='tight', dpi=400, transparent=True)
plt.show()
#%%
sns.boxplot(data=exclusive_df, linewidth=2)
plt.savefig("exclusive_ML_boxplot.png",bbox_inches='tight', dpi=400, transparent=True)
plt.show()

#%%

from scipy.stats import ranksums
ranksums(overlapped_df['G'], overlapped_df['MF'])


