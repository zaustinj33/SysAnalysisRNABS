import os
import glob
import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


# %%
# methods
# Read C_counts from fasta reference
def countChrs_from_file(input_idxstats):
    """
    Calculates the % of reads mapped to the mitochondria.

    Args:
        SAMtools idxstats of meRanGh files.
    Return:
        Raw counts and percentage of reads mapped to the mitochondria.
    """
    sample_dict = {}
    length_list = {}
    mapped_read_list = {}
    chrom_list = []
    with open(input_idxstats, "r") as input:
        input_stats = csv.reader(input, delimiter="\t")
        for line in input_stats:
            if 'random' not in line[0]:
                # chrom_list.append(line[0])
                length_list[line[0]] = int(line[1])
                mapped_read_list[line[0]] = int(line[2])

    # sample_dict['chrom'] = chrom_list
    sample_dict['length'] = length_list  # length, mapped reads
    sample_dict['mapped_reads'] = mapped_read_list

    return sample_dict


# %%

plot_dict = {}
norm_compare = {}
raw_compare = {}
raw_all = {}
ERCC_compare = {}
# Map_mapStats.txt for global map stats
# chrStats.txt for Cutoff analysis
for file in glob.glob("**/*chrStats.txt", recursive=True):  # remove '**/' if files are local to cwd
    name = os.path.basename(file).split("Map_mapStats.txt")[0]
    print(file)
    plot_dict[name] = countChrs_from_file(file)
    plot_df = pd.DataFrame.from_dict(plot_dict[name])

    # plot bargraph of total reads and normalized to chromosome length
    plot_df['normalized'] = plot_df['mapped_reads'] / plot_df['length']
    ERCC_df = plot_df.filter(like='ERCC', axis=0)
    ERCC_compare[name] = (ERCC_df['mapped_reads'].sum())
    plot_df = plot_df.filter(like='chr', axis=0)
    #print(plot_df['normalized']['chrM'])
    norm_compare[name] = plot_df['normalized']['chrM'] / sum(plot_df['mapped_reads'])+sum(ERCC_df['mapped_reads'])  # save for collective barplot
    raw_compare[name] = plot_df['mapped_reads']['chrM'] / sum(plot_df['mapped_reads'])+sum(ERCC_df['mapped_reads'])
    raw_all[name] = [plot_df['mapped_reads']['chrM'], sum(plot_df['mapped_reads'])-plot_df['mapped_reads']['chrM']]
    # unquote for individual plots

    # reorder xticks
    reorder_list = []
    for i in range(1,20):
        reorder_list.append("chr"+str(i))
    reorder_list = reorder_list + ["chrM", "chrX"]

    plot_df = plot_df.reindex(reorder_list)
    plot_df['perc_coverage_norm'] = plot_df['normalized'] / sum(plot_df['normalized'])
    plot_df['perc_coverage_raw'] = plot_df['mapped_reads'] / sum(plot_df['mapped_reads'])
    #print(plot_df)

    # add break to normalized graph
    #avg_cov_fig, avg_rate_ax1 = plt.subplots()
    #avg_rate_ax1.set(yscale="log")
    #sns.barplot(ax=avg_rate_ax1, x=plot_df.index, y='perc_coverage_norm', data=plot_df, color='grey')
    #plt.xticks(rotation=90)
    #plt.show()
    #plt.savefig(name+'normalized_mapped_perc.png', bbox_inches='tight', dpi=400, transparent=True)

    #avg_cov_fig2, avg_rate_ax2 = plt.subplots()
    #avg_rate_ax2.set(yscale="linear")
    #sns.barplot(ax=avg_rate_ax2, x=plot_df.index, y='perc_coverage_raw', data=plot_df, color='grey')
    #plt.xticks(rotation=90)
    #plt.savefig(name+'rawcounts_mapped_perc.png', bbox_inches='tight', dpi=400, transparent=True)

# %%
def format_df_MTanalysis(df):
    df = pd.DataFrame(df.items())
    df['map_approach'] = df[0].str.contains("TX")
    df['map_approach'] = np.where(df['map_approach'] == True, "Trans", "Genome")
    df['sample'] = df[0].str.replace(r'817.*$|_.*$', '')
    df = df.sort_values(0)
    return df


def format_df_CutoffAnalysis(df):
    df = pd.DataFrame(df.items())
    df['map_approach'] = df[0].str.contains("Ccutoff")
    df['map_approach'] = np.where(df['map_approach'] == True, "Cutoff_3C", "All_reads")
    df['sample'] = df[0].str.replace(r'817.*$|_.*$', '')
    df = df.sort_values(0)
    return df

# %%

# All samples' MT raw coverage
MT_df = format_df_CutoffAnalysis(raw_compare)
MT_df = MT_df.sort_values('map_approach')
#%%
order_x = ['G5', 'G1', 'G2','G3', 'G4', 'MFRNAseq', 'MF', 'pMFRNAseq', 'pMF', 'SRR']
order_MT = ['G1', 'G2','G3', 'G4']

p = sns.barplot(x='sample', y=1, data=MT_df, hue='map_approach',edgecolor="black", linewidth=2,
                palette=['purple', 'red'], order=order_MT)
sns.despine()
sns.set(rc={'figure.figsize':(5,5)})
sns.set_style("ticks")
plt.xticks(rotation=90)
plt.ylim(0,0.03)
p.set(yticks=np.arange(0, 0.031, 0.01))
#plt.legend([],[], frameon=False)

plt.savefig('rawcounts_all_perc_MTreads_Ccutoff.png', bbox_inches='tight', dpi=400, transparent=True)
plt.show()

#%%

# All samples' MT normalized coverage
norm_df = format_df(norm_compare)

f_norm, ax_norm = plt.subplots()
sns.barplot(x='sample', y=1, data=norm_df, hue='map_approach', ax=ax_norm)
plt.xticks(rotation=90)
plt.show()
plt.savefig('rawcounts_all_perc_MTmapped.png', bbox_inches='tight', dpi=400, transparent=True)
