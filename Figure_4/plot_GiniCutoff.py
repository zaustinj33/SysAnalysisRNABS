import os
import json
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

#%%
from scipy.spatial.transform import rotation

with open(os.getcwd()+'\\Gall_Ginivalues.json', 'r') as data:
    gini_counts = json.loads(data.read())

df = pd.DataFrame.from_dict(gini_counts, orient='index')
#df = df.replace(0, np.nan)
#df = df.replace(np.nan, 0)
#df = df.replace(0, np.nan)


mask = np.zeros_like(df)

#%%
order = ['G1','G2','G3', 'G4','MF_rep1', 'MF_rep2','pMF_rep1', 'pMF_rep2','SRR8170377', 'SRR8170378', 'SRR8170379', 'SRR8170380']

cols = df.columns.tolist()
#cols = cols[-1:] + cols[:-1]
df = df[cols]

sns.heatmap(df, annot=False, cmap="YlGnBu_r", linewidth=0.5,vmin=0, vmax=0.8)#, order=order)#, center=df.loc['MF_rep1', '2'])
plt.savefig("GiniCoef_heatmap.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()

plt.close()

#%%

sns.boxplot(data=df.T, linewidth=2)#, center=df.loc['MF_rep1', '2'])
#plt.savefig("GiniCoef_lineplot.png", bbox_inches='tight', dpi=400, transparent=True)
plt.xticks(rotation=90)
plt.show()

plt.close()
