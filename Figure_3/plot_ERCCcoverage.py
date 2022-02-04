import glob, os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

#%%

def read_ERCCmix():
    ERCC_file = glob.glob("**/*ERCCmix_conc.txt", recursive=True)[0]
    ERCC_conc = pd.read_csv(ERCC_file, sep="\t", usecols=[1,3], header=0,
                            names=['ERCC ID', 'concentration'])
    ERCC_conc = ERCC_conc.reindex(columns=['ERCC ID',  'concentration'])

    return ERCC_conc

#%%

def create_ERCC_cov():
    """
    Function to find SAMtools idxstats files of samples and record read coverage of ERCCs
    Return:
        ERCC counts for all relevant files with _mapStats
    """
    ERCC_list = []
    name_list = []
    df_headers = ['ERCC ID', 'length', 'coverage', 'unmapped_coverage']
    for file in glob.glob("**/*_mapStats.txt", recursive=True):
        print(file)  # name output file
        name = os.path.basename(file).split("_mapStats.txt")[0]
        print(name)
        df_temp = pd.read_csv(file, sep='\t', names = df_headers)
        df_temp = df_temp[df_temp['ERCC ID'].str.contains("ERCC")]
        name_list.append(name)
        ERCC_list.append(df_temp)

    result = pd.concat(ERCC_list, keys=name_list)
    return result


#%%
mix = read_ERCCmix()
coverage_df = create_ERCC_cov()
mix_sorted = mix.sort_values(by='concentration', ascending=True)
mix_sorted_name = mix_sorted['ERCC ID'].to_list()

coverage_df.index.names = ['name', 'anno']
# melt dataframe for better plotting
coverage_df = coverage_df.reset_index()
melt_df = pd.merge(coverage_df, mix, "left")
melt_df['concentration'] = melt_df['concentration'].astype(str)

melt_df['group'] = melt_df['name'].str.replace(r'817.*$|_.*$', '')  # add levels
melt_df['dup'] = np.where(melt_df['name'].str.contains('dedup'), "Dedup","Non-dedup")

# sort by coverage in mix
melt_df['ERCC ID'] = pd.Categorical(melt_df['ERCC ID'], mix_sorted_name)
melt_df = melt_df.sort_values('ERCC ID')
#melt_df = melt_df.replace(to_replace='External', value='Huang')

# or sort by length
melt_df['length'] = pd.Categorical(melt_df['length'])
melt_df = melt_df.sort_values('length')


#%%

order_hue = ['G5', 'G1', 'G2', 'G3', 'G4', 'SRR']
order_style = ['Non-dedup', 'Dedup']
#sns.stripplot(x='concentration', y='coverage', data=melt_df, hue='group', palette='bright',
 #             dodge=True, jitter=0.05, size=4, hue_order=order_hue)
sns.scatterplot(x='concentration', y='coverage', data=melt_df, ci=None, err_style="bars",
             hue='group', palette='bright', hue_order=order_hue, style='dup', style_order=order_style)
plt.yscale('log')
plt.xticks(rotation=90)
plt.legend([],[], frameon=False)

plt.savefig("Dotplot_ERCC_coverage.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()
plt.close()

