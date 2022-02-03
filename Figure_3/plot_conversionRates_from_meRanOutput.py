import glob, os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats

#%%

def conversion_rates():
    conv_names = []
    GN_conv = []
    TX_conv = []
    print()
    with open(glob.glob('**/conversion_rates.txt')[0], 'r') as file:
        input = [line.strip('/Total C to T conversion rate estimated:\t\n') for line in file]

    conv_names = input[0::3]
    conv_names_clean = [x[:len(x)-25] for x in conv_names]
    #conv_names_clean = [x.strip('/') for x in conv_names_clean]
    conv_names_GN = [x+"_GN" for x in conv_names_clean]
    conv_names_TX = [x+"_TX" for x in conv_names_clean]
    GN_conv = input[1::3]
    GN_conv_float = [float(x) for x in GN_conv]
    TX_conv = input[2::3]
    TX_conv_float = [float(y) for y in TX_conv]

    GN_df = pd.DataFrame(zip(conv_names_GN, GN_conv_float),
                           columns=['name', 'conv'])

    TX_df = pd.DataFrame(zip(conv_names_TX, TX_conv_float),
                         columns=['name', 'conv'])

    conv_df = pd.concat([GN_df, TX_df])

    return conv_df


test_df = conversion_rates()
#%%
test_df['group'] = test_df['name'].str.replace(r'817.*$|_.*$', '')  # add levels
test_df['map'] = test_df['name'].str.replace(r'.*_', '')

#%%
stats.ttest_ind(test_df['conv'][((test_df['map'] == 'GN') & (test_df['group'] == 'MF'))],
                test_df['conv'][((test_df['map'] == 'TX') & (test_df['group'] == 'MF'))])

#%%

sns.set_style("ticks")
sns.stripplot(x='group', y='conv', data=test_df, hue='map', split=True,
            palette='bright')
sns.boxplot(x='group', y='conv', data=test_df, hue='map', linewidth=2,
            palette='bright', showfliers=False)
sns.despine(offset=10, trim=True)
plt.legend([],[], frameon=False)

#plt.savefig("Overall_conv_rates.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()

plt.close()
