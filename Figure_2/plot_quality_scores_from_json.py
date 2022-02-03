"""
Interpreting fastp's json files for easy plotting of fq summary stats.
"""
import os
import glob
import json
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap


# %%

def get_nucleotide_quality(json_file):
    """
    Grabs mean quality score information from .json file and converts to plottable dictionary.

    Args:
        json_file: a .json file from fastQC software output.
    Return:
        1) dictionary of median nucleotide quality scores with the order:
        {[Filtering step] : {[R1 or R2] : {[Nucleotide] : [list of quality scores]}}}
    """
    test_json = json.loads(json_file.read())
    nuc_scores = {}
    if "length" in file:
        print("raw")
        nuc_scores['Total'] = {'R1_nuc_scores': test_json["read1_before_filtering"]["quality_curves"],
                               'R2_nuc_scores': test_json["read2_before_filtering"]["quality_curves"]}
        nuc_scores['Adapter'] = {'R1_nuc_scores': test_json["read1_after_filtering"]["quality_curves"],
                                 'R2_nuc_scores': test_json["read2_after_filtering"]["quality_curves"]}
        print(list(nuc_scores.keys()))

    else:
        print("processed")
        nuc_scores['R1_nuc_scores'] = (test_json["read1_after_filtering"]["quality_curves"])
        nuc_scores['R2_nuc_scores'] = (test_json["read2_after_filtering"]["quality_curves"])  # get processed reads

    # weight = test_json['summary']['after_filtering']['total_reads']

    return nuc_scores

# For plotting individual scores from file
def hist_from_jsonScores_dict(nuc_quality_scores):
    """
    Converts list of nucleotide scores to histogram format.

    Args:
        nuc_quality_scores: list of Phred33 quality scores.
    Return:
        Binned quality scores in percentile format.
    """
    test_df = pd.DataFrame.from_dict(nuc_quality_scores).values.T
    # bin Phred values
    binned_df = np.digitize(test_df, bins=[0, 20, 30, 35, 42])[:, :4]
    binned_df = pd.DataFrame(binned_df)
    binned_df.columns = ['A', 'T', 'C', 'G']

    # melt DF for easy plotting
    aggre = pd.melt(binned_df, value_vars=['A', 'T', 'C', 'G'], ignore_index=False)
    counts = aggre.groupby(['value', 'variable'])['value'].count().unstack()
    counts_perc = counts.div(counts.sum()).T  # convert to percentage
    # fill in columns with NA if missing
    for i in range(1, 5):
        if i not in counts_perc:
            counts_perc[i] = np.nan
    # sort columns by numeric
    counts_perc = counts_perc.reindex(columns=[1, 2, 3, 4])
    counts_perc = counts_perc.reindex(['A', 'T', 'G', 'C'])
    # print(counts_perc)
    return counts_perc


# Plot histogram of quality scores for R1 R2 separately
def plot_QualityHist(data_input):
    """
    Plots quality scores in the form of a stacked histogram

    Args:
        data_input: hist_from_jsonScores_dict() output; binned percentile dataframe
    Return:
        'qualityScore_stack.png'
    """
    global R1_scores, R2_scores
    for read_side, scores in data_input.items():
        read_name = read_side
        # print(read_name)

        hist_df = []
        for nuc, nuc_scores in scores.items():
            # print(nuc)
            hist_df.append(nuc_scores)  # merge nucleotide lists
            # together regardless of position
        # assign scores to read, plot paired plot
        if 'R1' in read_name:
            R1_scores = hist_from_jsonScores_dict(hist_df)
        elif 'R2' in read_name:
            R2_scores = hist_from_jsonScores_dict(hist_df)
        else:
            print("Read score not valid")

    custom_palette = ['#9a0007', '#d32f2f', '#ff6659', '#ef9a9a']
    color_dict = sns.color_palette(custom_palette)
    cmap = LinearSegmentedColormap.from_list("my_colormap", color_dict)

    bar_fig, (bar_ax1, bar_ax2) = plt.subplots(1, 2, sharex='all')
    R1_scores.plot(kind='bar', stacked=True, rot=1, ax=bar_ax1, colormap=cmap,
                   legend=None, width=1.0, edgecolor='black')
    R2_scores.plot(kind='bar', stacked=True, rot=1, ax=bar_ax2, colormap=cmap,
                   legend=None, width=1.0, edgecolor='black')
    bar_ax1.spines.top.set_visible(False)
    bar_ax2.spines.top.set_visible(False)
    bar_ax1.spines.right.set_visible(False)
    bar_ax2.spines.right.set_visible(False)
    bar_ax1.set_title(name)

    plt.savefig(os.getcwd() + "\\" + name + "_qualityScore_stack.png", bbox_inches='tight', dpi=400,
                transparent=True)
    plt.show()
    plt.close()
    # return hist_df


# %%
# main
"""
Iteration loop to aggregate data of all steps in fq filtering pipeline.

.json files should contain the suffixes:
'length'
'Quality filter'
'Window trimmed'
'm-bias trimmed'

    Args:
        all_BS_list.txt: names of input samples in new-line seperated file, eg:
        MF_batch1_rep1
        MF_batch1_rep2
        ...
    Return:
        complete_dict: a multiindex dataframe containing mean quality scores of all samples seperated by 
        step in the pipeline
"""

# Initialize dataframe from input list.txt
# Generate dictionary scores from files based on names in input.txt list
complete_dict = {}
# input_list = sys.argv[1]
with open("all_BS_list.txt", 'r') as data_list:
    for name in data_list:
        name = name.strip()
        # print(name)

        for file in glob.glob("**/" + name + "*.json", recursive=True):  # load json files by name
            plot_dict = {}
            try:
                with open(file, 'r') as j:
                    # By position, nucleotide mean quality score
                    full_name = os.path.basename(file).strip('.json')
                    if "length" in file:
                        temp_dict = get_nucleotide_quality(j)
                        plot_dict['Raw reads'] = temp_dict['Total']
                        complete_dict[full_name + "_Total"] = plot_dict

                        plot_dict = {'Adapter trimmed': temp_dict['Adapter']}  # reset dict for consistent dict indexing
                        complete_dict[full_name + "_Adapter"] = plot_dict
                    elif "qual" in file:
                        plot_dict['Quality filter'] = get_nucleotide_quality(j)  # or quality failed
                        complete_dict[full_name] = plot_dict
                    elif "Window" in file:
                        plot_dict['Window trimmed'] = get_nucleotide_quality(j)  # or window failed
                        complete_dict[full_name] = plot_dict
                    elif "4bpTrim" in file:
                        plot_dict['m-Bias trimmed'] = get_nucleotide_quality(j)  # or mbias trimmed failed
                        complete_dict[full_name] = plot_dict
                    else:
                        print("unknown source of file")
                        # plot_dict = {'nuc_quality_scores': get_nucleotide_quality(j)}
                    print(file)

            except EnvironmentError:
                print(file + " does not exist, need to create json file")

# %%

# Iteration loop to plot quality scores stacked bar for each sample file
for key, value in complete_dict.items():
    #   plot_quality_scores(complete_dict[key]['nuc_quality_scores'])
    full_name = key
    print(key)
    plot_QualityHist(complete_dict[key]['nuc_quality_scores'])
# %%

"""
Combine quality scores from clean reads from all libraries as Multi-index df
for comparing quality scores between files, between samples

Discerns between 'C' and 'G' bases on R1 and R2 as p-m5C sites
    
    Args:
        complete_dict: main() output; multiindex dataframe
    Return:
        test_complete: dictionary used to plot p-m5C quality vs other nucleotides as boxplot
"""

test_complete = {}
for sample, innerDict in complete_dict.items():
    # print(sample)
    for filter_type, innerValues in innerDict.items():
        # print(filter_type)
        for Read_side, values in innerValues.items():
            # print(Read_side)
            if "G" not in sample and Read_side == 'R1_nuc_scores':  # Swap read orientation to match MT libraries
                # have R1 = C, R2 = G
                Read_side = 'R2_nuc_scores'
            elif "G" not in sample and Read_side == 'R2_nuc_scores':
                Read_side = 'R1_nuc_scores'
            for nucleotide, scores in values.items():
                # fill in missing scores values if shorter than 150
                while len(scores) < 150:
                    scores.append(None)
                test_complete[(sample, Read_side, filter_type, nucleotide)] = scores


# %%

# Format test_complete
test_complete_df = pd.DataFrame.from_dict(test_complete).T
complete_melt = test_complete_df.stack(level=0).reset_index(level=4, drop=False).reset_index()
complete_melt['level_5'] = complete_melt['level_0'].str.replace(r'_failed.*$', '')
# %%

# Print summary statistics for sanity check
print("mean R1 score C \n", complete_melt[(complete_melt['level_3'] == 'C') &
                                          (complete_melt['level_2'] == 'm-Bias trimmed') &
                                          (complete_melt['level_1'] == 'R1_nuc_scores')].mean())
print("mean R1 score non-C \n", complete_melt[(complete_melt['level_3'] != 'C') &
                                              (complete_melt['level_2'] == 'm-Bias trimmed') &
                                              (complete_melt['level_1'] == 'R1_nuc_scores')].mean())

print("mean R2 score G \n", complete_melt[(complete_melt['level_3'] == 'G') &
                                          (complete_melt['level_2'] == 'm-Bias trimmed') &
                                          (complete_melt['level_1'] == 'R2_nuc_scores')].mean())
print("mean R2 score non-G \n", complete_melt[(complete_melt['level_3'] != 'G') &
                                              (complete_melt['level_2'] == 'm-Bias trimmed') &
                                              (complete_melt['level_1'] == 'R2_nuc_scores')].mean())
# %%

# Plot by p-m5C quality, separate R1 and R2
# order boxplot axis
complete_melt['level_2'] = pd.Categorical(complete_melt['level_2'],
                        ['Raw reads', 'Adapter trimmed',  'Window trimmed'], ordered=True)
# 'Quality filter', 'm-Bias trimmed'

def plot_m5C_boxplot(read, nucleotide, palette, axis):
    """
    Plots quality scores of nucleotides [A,T,C,G] as boxplot

    Args:
        read: Read 1 or Read 2
        nucleotide: 'A', 'T', 'C', or 'G'
        palette: Seaborn color palette
    Return:
        'qualityScore_stack.png'
    """
    sns.boxplot(x='level_5', y=0, hue='level_2', data=complete_melt.loc[
        ((complete_melt["level_3"] == nucleotide) &
         (complete_melt["level_1"] == read))],
                ax=axis, showfliers=False, palette=palette)
    axis.legend().remove()
    axis.set(xlabel=None)
    axis.set(ylabel=None)
    axis.set(ylim=[10,42])

# plots plot_m5C_boxplot() as multi-axis
fig, ((ax_C, ax_G), (ax_mean1, ax_mean2)) = plt.subplots(2, 2, sharex='all', sharey='all')
ax_C = plot_m5C_boxplot('R1_nuc_scores', 'C', 'OrRd', ax_C)
ax_G = plot_m5C_boxplot('R2_nuc_scores', 'G', 'OrRd', ax_G)
ax_mean1 = plot_m5C_boxplot('R1_nuc_scores', 'mean', 'OrRd', ax_mean1)
ax_mean2 = plot_m5C_boxplot('R2_nuc_scores', 'mean', 'OrRd', ax_mean2)

fig.set_size_inches(25,15)
plt.ylim(10, 42)
plt.xlabel(None)
plt.xticks(rotation=90)
plt.show()

#plt.savefig('QualityHistAll_R1R2_separate.png', bbox_inches='tight', dpi=400, transparent=True)
plt.close()
#%%

# Plot by p-m5C quality, combine R1 (C) and R2 (G) scores

def plot_m5C_boxplot_combine(nucleotide1, nucleotide2, palette, axis):
    """
    Plots quality scores comparing two nucleotides [A,T,C,G] as boxplot. Merges R1 and R2 information.

    Args:
        nucleotide1, nucleotide2: 'A', 'T', 'C', or 'G'
        palette: Seaborn color palette
        axis: seaborn axis
    Return:
        'qualityScore_stack.png'
    """
    sns.boxplot(x='level_5', y=0, hue='level_2', data=complete_melt.loc[
        ((complete_melt["level_3"] == nucleotide1) & (complete_melt["level_1"] == 'R1_nuc_scores')) |
        ((complete_melt["level_3"] == nucleotide2) & (complete_melt["level_1"] == 'R2_nuc_scores'))],
                ax=axis, showfliers=False, palette=palette)
    #axis.legend().remove()
    axis.set(xlabel=None)
    axis.set(ylabel=None)
    axis.set(ylim=[10,42])


# %%

fig_combine, (ax_pm5C, ax_mean_combine) = plt.subplots(2, sharex='all')
ax_pm5C = plot_m5C_boxplot_combine('C', 'G', 'OrRd', ax_pm5C)
ax_mean_combine = plot_m5C_boxplot_combine('mean', 'mean', 'OrRd', ax_mean_combine)

fig_combine.set_size_inches(20,10)
plt.ylim(10, 42)
plt.xlabel(None)
plt.xticks(rotation=90)
plt.show()

#plt.savefig('QualityHistAll_R1R2_combine.png', bbox_inches='tight', dpi=400, transparent=True)
plt.close()
