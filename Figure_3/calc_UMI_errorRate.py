import csv
import os
from collections import Counter
# %%
# Create C_counts from sam file
def countCs_from_SAM(sam_file):
    """
    Records the p-m5C coverage and position-specific conversion rate of a SAM file.

    Args:
        sam_file: A small SAM file containing ERCC-mapped reads only
    Return:
        Two dictionaries of nucleotide positions, one recording C-counts, the other C position coverage
    """
    read_dump = []
    test = sam_file
    UG_dict_coverage = {}
    UG_dict_counts = {}
    with open(test, 'r') as read_sam_file, \
            open(os.getcwd() + "\\" + test + "_withCTs.sam", "a") as write_sam_file, \
            open(os.getcwd() + "\\" + test + "_UMIgroups_allC.txt", "a") as write_UMIgroupFile:
        write_sam_file.truncate(0)  # delete old contents
        write_UMIgroupFile.truncate(0)

        csv_writer = csv.writer(write_sam_file, delimiter='\t', lineterminator='\n')

        # 1) find UMI groups with "C" (or other nucleotide of interest)
        for line in csv.reader(read_sam_file, delimiter='\t'):
            if len(line) == 20:  # check line has UMI group tag
                C_count_dict = create_counter_dict()

                ## Get UMI group ID
                try:
                    UG_group = line[18].strip('UG:i:')  # UMI ID

                    ## Update coverage of UMI group
                    if UG_group in UG_dict_coverage:
                        UG_dict_coverage[UG_group] += 1
                    else:
                        UG_dict_coverage[UG_group] = 1
                    # read = line[0]
                    # chrom = line[2]
                    # absolute_position_in_reference = int(line[3])
                    # start = absolute_position_in_reference
                    # end = start + len(line[9])

                    # Update counts of UMI group by read position
                    base_count = 0
                    for base in line[9]:  # split sequence to bases, capture any Cs
                        base_count += 1
                        if base == 'T':  # Can change this to get counts for any nucleotide
                            # print("HIT")
                            C_count_dict[base_count] += 1

                            # Option write line to file
                            csv_writer.writerow(line)
                    # print(C_count_dict)

                    # either sum hits (coverage) of group to dict, or create new entry
                    if UG_group in UG_dict_counts:
                        UG_dict_counts[UG_group] = dict(Counter(UG_dict_counts[UG_group]) + Counter(C_count_dict))
                    else:
                        UG_dict_counts[UG_group] = C_count_dict

                except IndexError:  # Reads with UMI, but tag is mangled (unlikely)
                    UG_group = 'Group_not_found'
                    # print(UG_group)

            else:  # Reads with no UMI
                read_dump.append(line)
        read_sam_file.close()
        """  ## Not working
        # 2) Write unique UMIs with only concordant Cs to file (non-converted)
        concordant_UMIs = []
        for k, v in UG_dict_counts.items():  # for each UMI group's coverage
            if list(v.values()):
                print(list(v.values()))
                value = list(v.values())[0]
                if value > 1:
                    if all(value == pos for pos in v.values()):  # if all counts in UMI group match coverage
                        # print('hit: ' + k)
                        write_UMIgroupFile.write(k + '\n')  # then UMI is 100% concordant on that base
                        concordant_UMIs.append(k)

        # print(concordant_UMIs)
        """
    return UG_dict_counts, UG_dict_coverage


def create_counter_dict():
    count_dict = {}
    for i in range(0, 151):
        count_dict[i] = 0
    return count_dict

# %%

# record_dict = {UG_group : {1 : 1_counters, 2 : 2_counters, 3 : 3_counters....150: 150_counters}}
test_counts, test_coverage = countCs_from_SAM('test_data/G4_2D_meRanGh_genomeMap_dedupGroup.bam_ERCC.sam')

# %%
rate_list = []
test_rate = test_counts
for k, v in test_counts.items():
    for i, j in v.items():
        if j >= 1 and test_coverage[k] >= 10:
            print("Group:" + k + ", Position : " + str(i) + ", Counts: " + str(int(j)) + ", Coverage : "
                  + str(test_coverage[k]) + ", Rate: " + str((int(j) / test_coverage[k]))
                  )
            rate_list.append((int(j) / test_coverage[k]))
            # test_rate[k][i] = int(j)/int(test_coverage[k])

# %%
import pandas as pd

rate_df = pd.DataFrame(rate_list)
# %%

import seaborn as sns
import matplotlib.pyplot as plt

sns.stripplot(data=rate_df, y=0)
plt.show()

