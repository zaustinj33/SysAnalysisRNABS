import csv
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import glob, os
import pickle
import numpy as np

# %%


# Create C_counts from sam file
def UMI_matrix_from_SAM(sam_file):
    read_dump = []
    test = sam_file
    UG_name_list = []
    UG_matrix = {}  # {UMI_ID: [Series of nucleotides in read, ..., ...]}
    with open(test, 'r') as read_sam_file: #, open(os.getcwd()+"\\"+test+"_withCs.sam", "a") as write_sam_file:

        # 1) compile read sequence into list of series
        for line in csv.reader(read_sam_file, delimiter='\t'):  #
            try:
                UG_ID = line[18]  # check line has UMI group tag
                #print(UG_ID)

                if UG_ID not in UG_name_list:  # If UMI group doesn't exist, create it and the first entry
                    UG_name_list.append(UG_ID)
                    UG_matrix[UG_ID] = [pd.Series(list(line[9]))]
                else:
                    UG_matrix[UG_ID].append(pd.Series(list(line[9])))

            except IndexError:
                #print("No UMI group or bad file format")
                pass

        #clean_UGs = {k:v for k,v in UG_matrix.items() if len(v) < 2}

        # 2) combine series into matrix of nucleotides, then count proportion of each
        UG_dfs = {}
        for ID, entries in UG_matrix.items():
            temp_df = pd.concat(entries, axis=1).T
            #if temp_df.value_counts
            temp_count = temp_df.apply(pd.Series.value_counts, normalize=True)
            #out_count = temp_count.dropna(thresh=len(temp_count) -2,  axis=1)

            UG_dfs[ID] = temp_count

        # 3) merge into multiindex df, stack by nucleotides
        UG = pd.concat(UG_dfs.values(), keys=UG_dfs.keys())
        data_plot = UG.stack().unstack(level=1)

        data_plot = data_plot.drop('N', axis=1, errors='ignore')   # drop N column
        data_plot = data_plot[data_plot.isnull().sum(axis=1) < 3]   # drop groups with no parity (discordant case 2)
        data_plot = data_plot[data_plot['C'] > 0]  # only positions with C

    return data_plot

#%%
"""
test = UMI_matrix_from_SAM('Chapter_4/test_data/G3_1C_meRanGh_genomeMap_dedupGroup.bam_ERCC.sam')
#test['sum'] = test.drop('C', axis=1).sum(axis=1)

sns.set(style='ticks')

sns.catplot(data=test, kind='box')
#plt.yscale('log')
plt.xticks(rotation=90)
plt.show()
"""


#%%
plot_list = []
name_list = []
for file in glob.glob("**/*_meRanGh_genomeMapTrim_dedupGroup.bam_ERCC.sam", recursive=True):
    print(file)
    name = os.path.basename(file).split("_meRanGh_genomeMapTrim_dedupGroup.bam_ERCC.sam")[0]
    plot_list.append(UMI_matrix_from_SAM(file))
    name_list.append(name)

#%%

final_df = pd.concat(plot_list, keys=name_list)
stack_df = final_df.stack(list(range(final_df.columns.nlevels))).reset_index()
stack_df['group'] = stack_df['level_0'].str.replace(r'_.*$', '')

#%%
# save dataframe as table csv
cwd = os.getcwd()
path = cwd + "/UMIgrouped_stats.pkl"
final_df.to_pickle(path)
#%%
plot_df = pickle.load(open("UMIgrouped_stats.pkl", 'rb'))
#%%

sns.set(style='ticks')
p = sns.FacetGrid(stack_df, col='group')
p.map(sns.boxplot, "level_3", 0)
#plt.yscale('log')
plt.xticks(rotation=90)
plt.show()

#%%

# Add group to multiindex
final_df['group'] = final_df.index.get_level_values(0).str.replace('_.*$', '')  # add levels
final_df.set_index('group', append=True, inplace=True)

def create_stackplot_df(input_df, slice):
    stackplot_df = input_df[input_df.index.get_level_values(3) == slice]

    def entropy(row):
        input = row.dropna()
        entr = -np.sum(input * np.log2(input))
        return entr

    stackplot_df['entropy'] = stackplot_df.apply(entropy, axis=1)
    stackplot_df = stackplot_df.sort_values(['A','G','C','T','entropy'],
                                          ascending=['False','False','False','False','False'])
    return stackplot_df.drop(['entropy'], axis=1)
#%%
# stacked area percentage

plot_list = ['G1', 'G2', 'G3', 'G4','G5']

def plot_stack(input_group):

    final_df_G1 = create_stackplot_df(final_df, input_group)

    color_dict = {'A':'blue','T':'orange','G':'green','C':'red'}
    ax = final_df_G1.plot(kind='area', stacked=True, color=color_dict, linewidth=0)
    ax.set_ylabel('Percent (%)')
    ax.margins(0, 0) # Set margins to avoid "whitespace"
    plt.ylim(0,1)

    plt.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=False)

    #plt.savefig(input_group+"_UMIstackplot.png", bbox_inches='tight', dpi=400, transparent=True)
    plt.show()
    plt.close()

#%%

for input_group in plot_list:
    print(input_group)
    plot_stack(input_group)
