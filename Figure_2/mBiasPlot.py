"""
Converts Bismark m-bias information to plottable information seperated by read and sample.
"""
# import numpy as np
import glob
import os
from itertools import islice
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
#%%
def Plot_mBias(input_bismark_stats):
    """
    Plots methylation bias information from Bismark bismark_methylation_extractor

    Args:
        bismark_methylation_extractor "M-bias.txt" file
    Return:
        1) m-bias plot of R1 and R2 of sample file, saved to working directory
        2) dataframe of R1 R2 methylation stats for aggregate plotting
    """

    name = os.path.basename(input_bismark_stats).split(".M-bias.txt")[0]
    # Read in File, write lines to list
    df_list = []
    CpG_names = []
    col_names = ["position", "count methylated", "count unmethylated", "% methylation", "coverage"]
    with open(input_bismark_stats, 'r') as data_file:
        # data_reader = csv.reader(data_file, delimiter='\t')
        for line in data_file:
            if line.startswith("C", 0):
                CpG_names.append(line.strip('\n'))  # get df names
                next(data_file)  # skip repeating '===' line break
                df_list.append([''.join(islice(data_file, 152))])
            else:
                pass  # print("no CpG files found")

    dfs = []
    # Convert lists to data frames for CpG, CHG, CHH
    for items in df_list:
        df = str(items).split('\\n')
        df = pd.DataFrame([x.split('\\t') for x in df], columns=col_names).iloc[1:138, :]  # convert to df and remove
        df = df.apply(pd.to_numeric)
        # duplicate col names
        # print(df)
        dfs.append(df)
    # print(dfs)

    # zip lists to dictionary
    sample_dict = dict(zip(CpG_names, dfs))

    ## R1 v R2
    colors = ['purple', 'orange']
    read_dfs_dict = {"R1": pd.concat(dfs[1:3]).groupby(level=0).sum(),  # sum methylation counts to calc global me
                    "R2": pd.concat(dfs[4:6]).groupby(level=0).sum()}
    read_dfs_dict['R1']['% methylation'] = read_dfs_dict['R1']['count methylated']/read_dfs_dict['R1']['count unmethylated']
    read_dfs_dict['R2']['% methylation'] = read_dfs_dict['R2']['count methylated']/read_dfs_dict['R2']['count unmethylated']
    read_dfs_dict['R1']['position'] = read_dfs_dict['R1']['position']/2
    read_dfs_dict['R2']['position'] = read_dfs_dict['R2']['position']/2

    read_dfs_colors = dict(zip(["R1", "R2"], colors))

    mbias_fig, (mbias_ax_R1, mbias_ax_R2) = plt.subplots(2,1, sharex='all')
    #mbias_ax_R1.spines['left'].set_visible(False)
    #mbias_ax_R2.spines['left'].set_visible(False)
    #mbias_ax.spines['top'].set_visible(False)
    for key, value in read_dfs_dict.items():
        print(key)
        print("min %: ", min(read_dfs_dict[key]['% methylation']), '\n',
              "max %: ", max(read_dfs_dict[key]['% methylation']))

    color_R1 = read_dfs_colors['R1']
    color_R2 = read_dfs_colors['R2']
    mbias_ax_R1.plot(read_dfs_dict['R1']['% methylation'], color=color_R1)
    mbias_ax_R2.plot(read_dfs_dict['R2']['% methylation'], color=color_R2)
    #plt.legend(key)
    plt.colormaps()
    plt.xlabel("Position")

    mbias_ax_R1.set_ylim(0,1)
    mbias_ax_R2.set_ylim(0,2.5)
    mbias_ax_R1.yaxis.tick_right()
    mbias_ax_R2.yaxis.tick_right()
    mbias_ax_R1.set_yticks(np.arange(0, 1.1, 0.25))
    mbias_ax_R2.set_yticks(np.arange(0, 2.1, 0.5))
    mbias_ax_R1.set_ylim([-.1,1])
    mbias_ax_R2.set_ylim([-.1,2.1])
    #plt.show()
    plt.savefig(name +"_bothReads"+ '_mbias_trim.png', bbox_inches='tight', dpi=400, transparent=True)
    plt.close()


    return read_dfs_dict

# %%
# main
"""
Iteration loop over working directory to find processed mbias files and get stats
"""
merged_results_dict = {}
for file in glob.glob("**/G[3|4]*_1_20C_cutoffClean_bismark_bt2_pe.M-bias.txt", recursive=True):
    name = os.path.basename(file).split("_1_20C_cutoffClean_bismark_bt2_pe.M-bias.txt")[0]
    print(file)
    # gets stats
    results = Plot_mBias(file)
    merged_results_dict[name] = results
#%%
"""
Convert list of dictionaries of dictionaries to dataframe, then plot average methylation information
"""

# Initialize plot
mbias_fig, mbias_ax = plt.subplots()
mbias_ax.spines['right'].set_visible(False)
mbias_ax.spines['top'].set_visible(False)
colors = sns.color_palette(None, len(merged_results_dict))
#print(merged_results_dict)

# Create color dictionary based on number of samples
read_dfs_colors = {}
i = 0
for name, sample in merged_results_dict.items():
    read_dfs_colors[name] = colors[i]
    i += 1
#print(read_dfs_colors)

# Plot the (average?) % methylation at each position of read, for all samples
for key, value in merged_results_dict.items():
    print(key)
    color = read_dfs_colors[key]

    # plot average methylation, or R1 only methylation? Unhash for average methylation
    avg_meth = (value['R1']['% methylation'].astype(float) + value['R2']['% methylation'].astype(float))/2
    # mbias_ax.plot(value['R1']['% methylation'].astype(float), color=color)
    # mbias_ax.plot(value['R2']['% methylation'].astype(float), color=color)
    mbias_ax.plot(avg_meth, color = color)

# plt.legend(merged_results_dict.keys())
# plt.colormaps()
plt.xlabel("Position")
plt.ylabel("% Methylation")
plt.savefig('Test_mbias_R1.png', bbox_inches='tight', dpi=400, transparent=True)
plt.show()

