import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

#%%

def file_path_checker(name, suffix):
    print(name + " " + suffix)
    try:
        potential_paths = glob.glob("**/"+name+"_Genome10xCall*"+suffix+".txt", recursive=True)
        print(potential_paths)
        if len(potential_paths) < 1:  # No file exists
            print("no matches")
            return None
        elif len(potential_paths) > 1:  # multiple matches on wildcard, return first file found
            print("multiple matches found")
            return potential_paths[0]
        else:
            return potential_paths[0]

    except:
        print("no matches")
        return None


def line_count(input):
    try:
        df = pd.read_csv(input, sep='\t', low_memory=False)  # file

        # MT sites only
        #df = df[df["#SeqID"] == 'MT']

        line_counter = len(df.index) - 1  # minus header
        gene_counter = len(df.gene.value_counts()) - 1  # minus no feature genes
        cov_counter = np.sum(df['cov']) + np.sum(df['C_count'])

    except (TypeError, ValueError):  # if no file present
        print("input has a formatting error, or does not exist")
        return np.nan,np.nan,np.nan

    return line_counter, gene_counter, cov_counter


def divide_tuples(tup1, tup2):
    try:
        res = tuple(round((ele1 / ele2),6) for ele1, ele2 in zip(tup1, tup2))
    except RuntimeWarning:  # divide by 0
        res = (0,0,0)
    return res

def flatten_list(nest):
    input_list = []
    for item in nest:
        if isinstance(item, (list, tuple, set)):
            input_list.extend(flatten_list(item))
        else:
            input_list.append(item)
    return input_list

sample_list = ['G1','G2','G3','G4', #['G1', 'G1_1A', 'G1_2A', 'G2', 'G2_1B', 'G2_2B','G3', 'G3_1C', 'G3_2C','G4', 'G4_1D', 'G4_2D',]
               'MF_rep1', 'MF_rep2','pMF_rep1', 'pMF_rep2', 'SRR8170377', 'SRR8170378', 'SRR8170379', 'SRR8170380']

filter_list = ['10x','basicFilter', 'Ccutoff', 'Cutoff+basicFilter', 'signalNoise', 'fdrFilter', 'strucFilter'] #
feature_list = ['site_count', 'gene_count', 'coverage_count']
# build multiindex from data lists
index = []
for i in sample_list:
    for j in filter_list:
        for k in feature_list:
            index.append((i,j,k))
#%%
## Percent of 10x sites
data_list = []
for name in sample_list:
    #try:
        sites_10x = (line_count(file_path_checker(name, "_annotateALL")))  # sites after 10x only
        sites_BF = divide_tuples(line_count(file_path_checker(name, "_basicFilterALL")), sites_10x)
        sites_C = divide_tuples(line_count(file_path_checker(name, "Cutoff_annotate")), sites_10x)  # sites after C-cutoff
        sites_BFC = divide_tuples(line_count(file_path_checker(name, "Cutoff_basicFilter")), sites_10x)  # sites after basic filter
        sites_DF = divide_tuples(line_count(file_path_checker(name, "signalFilter")), sites_10x)  # after signal Noise filter
        sites_FDR = divide_tuples(line_count(file_path_checker(name, "fdrFilter")), sites_10x)  # after signal FDR filter
        sites_STR = divide_tuples(line_count(file_path_checker(name, "strucFilter")), sites_10x)  # after signal Noise filter

        data_list.append([(1,1,1), sites_C, sites_BF, sites_BFC, sites_DF, sites_FDR, sites_STR])  #

    #except:
     #   print("mystery problem")
      #  pass

#%%

## Raw numbers
data_list = []
for name in sample_list:
    #try:
    sites_10x = (line_count(file_path_checker(name, "annotateALL")))  # sites after 10x only
    sites_BF = line_count(file_path_checker(name, "_basicFilterALL"))
    sites_C = line_count(file_path_checker(name, "Cutoff_annotate"))  # sites after C-cutoff
    sites_BFC = line_count(file_path_checker(name, "Cutoff_basicFilter")) # sites after basic filter
    sites_DF = line_count(file_path_checker(name, "signalFilter"))  # after signal Noise filter
    sites_FDR = line_count(file_path_checker(name, "fdrFilter"))  # after signal FDR filter
    sites_STR = line_count(file_path_checker(name, "strucFilter")) # after signal Noise filter

    data_list.append([sites_10x,sites_BF, sites_C, sites_BFC, sites_DF, sites_FDR, sites_STR]) #

#except:
#   print("mystery problem")
#  pass
#%%
flat = pd.Series(flatten_list(data_list), index=index)
MI_index = pd.MultiIndex.from_tuples(index, names=('sample','filter','feature'))
plot_df = flat.reindex(MI_index)

#%%
#plot results of filtering

#1) unique sites
plot_sites = plot_df[:,:,'site_count'].reset_index()
plot_sites = plot_sites.rename({0:'value'}, axis=1)
plot_sites.loc[plot_sites['value'] <= 0, 'value'] = np.nan
plot_sites = plot_sites[plot_sites['filter'] != 'Ccutoff']

#sns.pointplot(data=plot_sites, x='filter', y=0, color='black', ci='sd')
sns.stripplot(data=plot_sites, x='filter', y='value', hue='sample', palette='icefire')
plt.legend('')
sns.lineplot(data=plot_sites, x='filter', y='value', hue='sample', palette='icefire', linewidth=2)
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
plt.yscale('log')
#plt.ylim(0,1.2)
plt.savefig("Perc_remaining_sites.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()
#%%
#2) unique genes
plot_sites = plot_df[:,:,'gene_count'].reset_index()
plot_sites = plot_sites.rename({0:'value'}, axis=1)
plot_sites.loc[plot_sites['value'] <= 0, 'value'] = np.nan
plot_sites = plot_sites[plot_sites['filter'] != 'basicFilter']
#sns.pointplot(data=plot_sites, x='filter', y=0, color='black', ci='sd')
sns.stripplot(data=plot_sites, x='filter', y='value', hue='sample', palette='icefire')
plt.legend('')
sns.lineplot(data=plot_sites, x='filter', y='value', hue='sample', palette='icefire', linewidth=2)
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
#plt.yscale('log')
plt.savefig("Perc_remaining_genes_MT.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()

#%%
#2) remaining read coverage
plot_sites = plot_df[:,:,'coverage_count'].reset_index()
plot_sites = plot_sites.rename({0:'value'}, axis=1)
plot_sites.loc[plot_sites['value'] <= 0, 'value'] = np.nan
plot_sites = plot_sites[plot_sites['filter'] != 'basicFilter']

#sns.pointplot(data=plot_sites, x='filter', y=0, color='black', ci='sd')
sns.stripplot(data=plot_sites, x='filter', y='value', hue='sample', palette='icefire')
plt.legend('')
sns.lineplot(data=plot_sites, x='filter', y='value', hue='sample', palette='icefire', linewidth=2)
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
plt.yscale('log')
plt.savefig("Perc_remaining_coverage_MT.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()
