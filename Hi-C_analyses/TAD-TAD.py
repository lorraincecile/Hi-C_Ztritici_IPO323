#%%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

#%%
name = "IPO323_WT"

# Load the .tsv file
tsv_df = pd.read_csv(f'/Users/iglavinchesk/WP1/WP1.2/data/TAD-TAD/{name}.filtered.tsv', sep='\t', header=0)
print(tsv_df)
# Load the .bed file
bed_df = pd.read_csv(f'/Users/iglavinchesk/WP1/WP1.2/data/TADs/{name}_bin3kb_norm_ICE_5_DEFAULT_domains.bed', sep='\t', header=None)

# Assuming the fourth column in the .bed file is the one you want to merge
bed_df.columns = ['chrom', 'start', 'end', 'ID','score','strand','thickStart','thickEnd','itemRgb']

# Merge the DataFrames on a common column, if applicable
# Here, we assume 'chrom' and 'start' columns are common for merging
merged_df = pd.merge(tsv_df, bed_df[['chrom', 'start', 'ID']].rename(columns={'ID': 'ID1'}), left_on=['chrom1', 'start1'], right_on=['chrom', 'start'], how='left') \
              .merge(bed_df[['chrom', 'start', 'ID']].rename(columns={'ID': 'ID2'}), left_on=['chrom2', 'start2'], right_on=['chrom', 'start'], how='left')
# Remove unnecessary columns
merged_df.drop(columns=['chrom_x', 'start_x', 'chrom_y', 'start_y'], inplace=True)

#%%
# Remove rows where ID1 and ID2 are the same
filtered_df = merged_df[merged_df["ID1"] != merged_df["ID2"]]

intra_df = filtered_df[filtered_df["chrom1"] == filtered_df["chrom2"]]
inter_df = filtered_df[filtered_df["chrom1"] != filtered_df["chrom2"]]

filtered_df.to_csv(f'WP1.2/data/TAD-TAD/{name}.TAD-TAD.tsv', sep='\t', header=False, index=False)
intra_df.to_csv(f'WP1.2/data/TAD-TAD/{name}.intraTAD-TAD.tsv', sep='\t', header=False, index=False)
inter_df.to_csv(f'WP1.2/data/TAD-TAD/{name}.interTAD-TAD.tsv', sep='\t', header=False, index=False)

# %%
print(len(intra_df))
print(len(inter_df))
# %%
# Load the HI-TADs file
hi_tads_df = pd.read_csv(f'/Users/iglavinchesk/WP1/WP1.2/data/HI-TADs/{name}_HI-TADs.bed', sep='\t', header=None)
hi_tads_df.columns = ['chrom', 'start', 'end', 'id']

# Label each row in filtered_df
def label_row(row, hi_tads_ids):
    id1_in_hi_tads = row['ID1'] in hi_tads_ids
    id2_in_hi_tads = row['ID2'] in hi_tads_ids
    if id1_in_hi_tads and id2_in_hi_tads:
        return 'HI-TAD x HI-TAD'
    elif id1_in_hi_tads or id2_in_hi_tads:
        return 'HI-TAD x R-TAD'
    else:
        return 'R-TAD x R-TAD'

hi_tads_ids = set(hi_tads_df['id'])
intra_df['label'] = intra_df.apply(label_row, axis=1, hi_tads_ids=hi_tads_ids)
intra_df['log2'] = np.log2(intra_df['odds_ratio'])
intra_df.to_csv(f'WP1.2/data/TAD-TAD/{name}.intraTAD-TAD.tsv', sep='\t', header=False, index=False)

inter_df['label'] = inter_df.apply(label_row, axis=1, hi_tads_ids=hi_tads_ids)
inter_df['log2'] = np.log2(inter_df['odds_ratio'])
inter_df.to_csv(f'WP1.2/data/TAD-TAD/{name}.interTAD-TAD.tsv', sep='\t', header=False, index=False)

filtered_df['label'] = filtered_df.apply(label_row, axis=1, hi_tads_ids=hi_tads_ids)
filtered_df['log2'] = np.log2(filtered_df['odds_ratio'])
filtered_df.to_csv(f'WP1.2/data/TAD-TAD/{name}.TAD-TAD.tsv', sep='\t', header=False, index=False)


# %%
df = intra_df

# Create the boxplot
plt.figure(figsize=(10, 6))
boxplot = sns.boxplot(x='label', y='log_ratio', data=df, patch_artist=True)

# Customize the box colors and hatches
for patch, label in zip(boxplot.artists, df['label'].unique()):
    if label == 'HI-TAD x HI-TAD':
        patch.set_facecolor('lightgray')
        patch.set_edgecolor('black')
    elif label == 'R-TAD x R-TAD':
        patch.set_facecolor('white')
        patch.set_edgecolor('black')
    elif label == 'HI-TAD x R-TAD':
        patch.set_facecolor('lightgray')
        patch.set_edgecolor('black')
        patch.set_hatch('//')

plt.title('Boxplot of log_ratio by Label')
plt.xlabel('Label')
plt.ylabel('log_ratio')

# Perform pairwise comparisons
labels = df['label'].unique()
comparisons = [(labels[i], labels[j]) for i in range(len(labels)) for j in range(i+1, len(labels))]
p_values = []
for (label1, label2) in comparisons:
    group1 = filtered_df[filtered_df['label'] == label1]['log_ratio']
    group2 = filtered_df[filtered_df['label'] == label2]['log_ratio']
    stat, p = mannwhitneyu(group1, group2)
    p_values.append((label1, label2, p))

# Annotate the plot with significance stars at different heights
def get_significance_label(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return 'ns'

y_max = filtered_df['log2'].max()
for i, (label1, label2, p) in enumerate(p_values):
    y, h, col = y_max + (i+1)*0.5, 0.2, 'k'
    x1, x2 = labels.tolist().index(label1), labels.tolist().index(label2)
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
    plt.text((x1+x2)*.5, y+h, get_significance_label(p), ha='center', va='bottom', color=col)

plt.show()
# %%
