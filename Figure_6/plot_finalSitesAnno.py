import glob
import json, os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams

#%%

with open(os.getcwd()+'\\Chapter_5\\AggregateMapReads.json', 'r') as data:
    map_counts = json.loads(data.read())

map_df = pd.DataFrame.from_dict(map_counts, orient='index')
map_df['group'] = map_df.index.str.replace(r'_.*$', '')
#%%
conditions = ["G","MF","pMF","SRR"]

site_count = {}
for set in conditions:
    for file in glob.glob("**/"+set+"_allUnion.csv", recursive=True):
        df = pd.read_csv(file, sep='\t', low_memory=False)

        site_count[set] = len(df)


site_df = pd.DataFrame(site_count.items())
map_mean = map_df.groupby('group').sum().reset_index()
site_norm = pd.DataFrame(site_df[1].div(map_mean[0])*100000)
site_norm['group'] = map_mean['group']
#%%

fig, (ax1, ax2) = plt.subplots(ncols=2)
sns.barplot(data=site_norm, x='group', y=0, ax=ax1, linewidth=2, edgecolor='black', palette='bright')
sns.swarmplot(data=map_df, x='group', y=0, size=5, ax=ax2, linewidth=2, color='black')
ax2.yaxis.set_ticks_position('right')
rcParams['figure.figsize'] = 10,4

plt.savefig("final_siteUnionCount_mapCounts.png",bbox_inches='tight', dpi=400, transparent=True)

plt.show()
plt.close()

