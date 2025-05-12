#%%
import pandas as pd
import pybedtools
import numpy as np

name = "IPO323+dim2"

# Load BED files into pandas
boundaries = pd.read_csv(f"/Users/iglavinchesk/WP1/WP1.2/data/TADs/{name}_bin3kb_norm_ICE_5_DEFAULT_boundaries.bed", sep="\t", header=None, names=["chrom", "start", "end", "boundary_id","boundary_score","x1"])
boundaries = boundaries[['chrom', 'start', 'end', 'boundary_id', 'boundary_score']]
domains = pd.read_csv(f"/Users/iglavinchesk/WP1/WP1.2/data/TADs/{name}_bin3kb_norm_ICE_5_DEFAULT_domains.bed", sep="\t", header=None, names=["chrom", "start", "end", "id","x1","x2","x3","x4", "x5"])
domains = domains[["chrom", "start", "end", "id"]]
scores = pd.read_csv(f"/Users/iglavinchesk/WP1/WP1.2/data/TADs/{name}_bin3kb_norm_ICE_5_DEFAULT_score.bedgraph", sep="\t", header=None, names=["chrom", "start", "end", "insulation_score"])


#%%
# Convert to BedTools objects
boundaries_bed = pybedtools.BedTool.from_dataframe(boundaries[["chrom", "start", "end", "boundary_id", "boundary_score"]])
domains_bed = pybedtools.BedTool.from_dataframe(domains[['chrom', 'start', 'end','id']])
scores_bed = pybedtools.BedTool.from_dataframe(scores)

#%%
# Find closest boundary for each TAD
closest_boundaries = boundaries_bed.closest(domains_bed, d=True).to_dataframe()
closest_boundaries.columns = ["chrom", "start", "end", "boundary_id", "boundary_score", "tad_chrom", "tad_start", "tad_end", "id","distance"]

# Compute mean insulation score for each TAD
tad_scores = pybedtools.BedTool.from_dataframe(domains).map(scores_bed, c=4, o="max").to_dataframe()
print(tad_scores.head())

tad_scores.columns = ["chrom", "start", "end", "tad_id", "avg_insulation_score"]

# Merge insulation scores and boundaries
merged = tad_scores.merge(closest_boundaries, left_on=["chrom", "start", "end"], right_on=["tad_chrom", "tad_start", "tad_end"])
print(merged)
merged = merged[["tad_chrom", "tad_start", "tad_end", "tad_id","avg_insulation_score", "boundary_score"]]

# Step 1: Calculate percentiles
insulation_threshold = np.average(merged["avg_insulation_score"])  # Top 10%
boundary_threshold = np.average(merged["boundary_score"])  # Bottom 10%
merged["insulation_diff"] = merged["avg_insulation_score"] - merged["boundary_score"]
difference_threshold = np.median(merged["insulation_diff"])  # Bottom 10%
difference_threshold2 = merged["insulation_diff"].quantile(0.80)

# Step 2: Filter TADs based on thresholds
filtered_tads = merged[
    (merged["avg_insulation_score"] >= insulation_threshold) &  
    (merged["boundary_score"] <= boundary_threshold) &
    (merged["insulation_diff"] >= difference_threshold) &
    (merged["insulation_diff"] >= difference_threshold2) 
]

# Step 3: Drop duplicates to keep unique TADs
filtered_tads = filtered_tads.drop_duplicates(subset=["tad_id"])
#%%
filtered_tads["type"] = "HI-TAD"
# Extract relevant columns for BED file (chrom, start, end, id)
bed_columns = filtered_tads[["tad_chrom", "tad_start", "tad_end", "tad_id", "type"]]  # Replace "chrom_x", "start_x", "end_x", "id" if necessary
# Write the DataFrame to a BED file
bed_columns.to_csv(f"/Users/iglavinchesk/WP1/WP1.2/data/HI-TADs/{name}_HI-TADs.bed", sep="\t", header=False, index=False)

bed_columns_non_filtered = merged[~merged["tad_id"].isin(filtered_tads["tad_id"])]  # TADs that did not pass
# Extract the same relevant columns for TADs that did not pass
bed_columns_non_filtered["type"] = "R-TAD"
bed_columns_non_filtered = bed_columns_non_filtered[["tad_chrom", "tad_start", "tad_end", "tad_id", "type"]]
bed_columns_non_filtered.to_csv(f"/Users/iglavinchesk/WP1/WP1.2/data/HI-TADs/{name}_R-TADs.bed", sep="\t", header=False, index=False)

# %%
print(filtered_tads[filtered_tads["tad_chrom"] == 4])
# %%
