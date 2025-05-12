#%%
import cooler
import numpy as np
import pandas as pd

name = 'IPO323+dim2'

# Load the .cool file
cool_file = f'/Users/iglavinchesk/WP1/WP1.2/data/bin10kb/{name}_bin10kb_norm_ICE_OE.cool'
c = cooler.Cooler(cool_file)

# Extract the matrix
matrix = c.matrix(balance=False)[:]

resolution = 10000
res_string = str(resolution)
chromosomes = c.chromnames
chr_starts = []
chr_ends = []
for chrom in chromosomes:
	chr_starts.append(c.extent(chrom)[0])
	chr_ends.append(c.extent(chrom)[1])

#%%
matrix = c.matrix(balance=False)[:]
chrom_array = np.array(matrix)
array_x_dim, array_y_dim = chrom_array.shape

# %%

def get_chrom_and_coordinates(index, chr_starts, chr_ends, resolution, chromosomes):
	for k in range(len(chr_starts)):
		if chr_starts[k] <= index < chr_ends[k]:
			chrom_string = chromosomes[k]
			start = (index - chr_starts[k]) * resolution
			end = start + resolution
			return chrom_string, start, end
	return None, None, None

def filter_value(value):
    """
    Filters out NaN, Inf, and zero values.
    
    Parameters:
    value (float): The value to be filtered.
    
    Returns:
    bool: True if the value should be skipped, False otherwise.
    """
    if np.isnan(value) or np.isinf(value) or value == 0.00:
        return True
    return False

    """
    Processes the chrom_array to log valid values and outputs the coordinates of interactions.

    Parameters:
	chrom_array (numpy array): The interaction matrix.
    chr_starts (list): List of chromosome start positions.
    chr_ends (list): List of chromosome end positions.
    resolution (int): The resolution for coordinates.
    chromosomes (list): Chromosome labels.

    Returns:
    pd.DataFrame: A DataFrame containing the results.
    """

def merge_adjacent_bins(df):
    """
    Merges adjacent bins in the DataFrame.
    
    Parameters:
    df (pd.DataFrame): The DataFrame containing the bins.
    
    Returns:
    pd.DataFrame: A DataFrame with merged adjacent bins.
    """
    merged_results = []
    for chrom1 in df['Chromosome1'].unique():
        chrom_df = df[df['Chromosome1'] == chrom1].sort_values(by='Start1')
        current_start1, current_end1 = chrom_df.iloc[0]['Start1'], chrom_df.iloc[0]['End1']
        current_chrom2, current_start2, current_end2 = chrom_df.iloc[0]['Chromosome2'], chrom_df.iloc[0]['Start2'], chrom_df.iloc[0]['End2']
        for _, row in chrom_df.iterrows():
            if row['Start1'] <= current_end1 and row['Chromosome2'] == current_chrom2 and row['Start2'] <= current_end2:
                current_end1 = max(current_end1, row['End1'])
                current_end2 = max(current_end2, row['End2'])
            else:
                merged_results.append([chrom1, current_start1, current_end1, current_chrom2, current_start2, current_end2])
                current_start1, current_end1 = row['Start1'], row['End1']
                current_chrom2, current_start2, current_end2 = row['Chromosome2'], row['Start2'], row['End2']
        merged_results.append([chrom1, current_start1, current_end1, current_chrom2, current_start2, current_end2])
    return pd.DataFrame(merged_results, columns=['Chromosome1', 'Start1', 'End1', 'Chromosome2', 'Start2', 'End2'])
#%%
results = []
x = 0
for i in range(array_y_dim):
    chrom_string1, start1, end1 = get_chrom_and_coordinates(i, chr_starts, chr_ends, resolution, chromosomes)
        
    for j in range(array_x_dim):
        chrom_string2, start2, end2 = get_chrom_and_coordinates(j, chr_starts, chr_ends, resolution, chromosomes)
        chrom1_int = int(chrom_string1)
        chrom2_int = int(chrom_string2)
        value = float(chrom_array[i, j])
            
        # Filter invalid values
        if filter_value(value):
            continue
        #print(chrom1_int, chrom2_int) 
        # Process valid values with log2 transform and threshold check
        if (value >= 3.5) and (chrom1_int < chrom2_int):
            x += 1
            print(x)
            # Convert chromosome strings to integers for comparison
            results.append([chrom_string1, start1, end1, chrom_string2, start2, end2])
    
    # Convert results to DataFrame
df = pd.DataFrame(results, columns=['Chromosome1', 'Start1', 'End1', 'Chromosome2', 'Start2', 'End2'])

print(df.head())
#%%
# Merge adjacent bins
merged_df = merge_adjacent_bins(df)
df = pd.DataFrame(results, columns=['Chromosome1', 'Start1', 'End1', 'Chromosome2', 'Start2', 'End2'])

# Save merged_df to a links file for Circos
df.to_csv(f'/Users/iglavinchesk/WP1/WP1.2/config/{name}_all_links.txt', sep='\t', index=False, header=False)
print(df)
# %%
# Define file path
cent_array = np.loadtxt('/Users/iglavinchesk/WP1/WP1.2/data/annotations/IPO323_centromeres3c.bed', dtype=int)

# Cache chromosome coordinates outside the loop
chrom_coords = [get_chrom_and_coordinates(i, chr_starts, chr_ends, resolution, chromosomes) for i in range(array_y_dim)]
x = 0 
results = []
# Loop through the array
for i in range(array_y_dim):
    chrom_string1, start1, end1 = chrom_coords[i]

    for j in range(array_x_dim):
        chrom_string2, start2, end2 = chrom_coords[j]
        
        # Fetch value and check if it is finite and non-zero
        value = float(chrom_array[i, j])
        if not np.isfinite(value) or value == 0.00:
            continue

        # Compare chromosomes and handle core-chromosome conditions
        if (value >= 3.5) & (chrom_string1 < chrom_string2):
            cent_start1, cent_end1 = map(int, cent_array[int(chrom_string1)-1, 1:3])  # Only get start & end
            cent_start2, cent_end2 = map(int, cent_array[int(chrom_string2)-1, 1:3])  # Only get start & end
            # Only compute range check once per loop iteration
            # Check if start1 and start2 are within 50kb of the centromere
            if (cent_start1 - 70000 <= start1 <= cent_start1 + 70000 or cent_end1 - 70000 <= start1 <= cent_end1 + 70000) and \
                (cent_start2 - 70000 <= start2 <= cent_start2 + 70000 or cent_end2 - 70000 <= start2 <= cent_end2 + 70000):
                results.append([chrom_string1, start1, end1, chrom_string2, start2, end2, "color=black_a5", "order1"])
                x += 1
                print(x)
# Otherwise, add all other interactions
            else:
                results.append([chrom_string1, start1, end1, chrom_string2, start2, end2, "color=color_13_a5", "order0"])
                x += 1

# %%
df = pd.DataFrame(results, columns=['Chromosome1', 'Start1', 'End1', 'Chromosome2', 'Start2', 'End2','Color', "Order"])
print(df.head())
df = df.sort_values(by="Order")
# Save merged_df to a links file for Circos
df.to_csv(f'/Users/iglavinchesk/WP1/WP1.2/config/{name}_CentromereColored_links.txt', sep='\t', index=False, header=False)

# %%
# Define file path

core = range(1,14)
accessory = range(14, 21)

# Cache chromosome coordinates outside the loop
chrom_coords = [get_chrom_and_coordinates(i, chr_starts, chr_ends, resolution, chromosomes) for i in range(array_y_dim)]
x = 0 
results = []
# Loop through the array
for i in range(array_y_dim):
    chrom_string1, start1, end1 = chrom_coords[i]

    for j in range(array_x_dim):
        chrom_string2, start2, end2 = chrom_coords[j]
        
        # Fetch value and check if it is finite and non-zero
        value = float(chrom_array[i, j])
        if not np.isfinite(value) or value == 0.00:
            continue

        # Compare chromosomes and handle core-chromosome conditions
        if (value >= 3.5) & (chrom_string1 < chrom_string2):
            # Only compute range check once per loop iteration
            # Check if start1 and start2 are within 50kb of the centromere
            if (int(chrom_string1) in core and int(chrom_string2) in accessory) or (int(chrom_string1) in accessory and int(chrom_string2) in core):
                results.append([chrom_string1, start1, end1, chrom_string2, start2, end2, "color=color_5_a5", "order1"])
                x += 1
                print(x)
# Otherwise, add all other interactions
            else:
                results.append([chrom_string1, start1, end1, chrom_string2, start2, end2, "color=color_13_a5", "order0"])
                x += 1

# %%
df = pd.DataFrame(results, columns=['Chromosome1', 'Start1', 'End1', 'Chromosome2', 'Start2', 'End2','Color', "Order"])
print(df.head())
df = df.sort_values(by="Order")
# Save merged_df to a links file for Circos
df.to_csv(f'/Users/iglavinchesk/WP1/WP1.2/config/{name}_CoreAccessory_links.txt', sep='\t', index=False, header=False)

# %%

# Define bin size (10 kb)
bin_size = 10000
k9_array = np.transpose(np.genfromtxt(
    #'/Users/iglavinchesk/WP1/WP1.2/data/ChIPseq/WT_H3K9me3_mean.10kb.clean.bed',
    '/Users/iglavinchesk/WP1/WP1.2/data/5mC/IPO323+dim2_5mC_w0.10kb.bed',
    usecols=3, 
    delimiter=None,  # Adjust if needed
    dtype=float, 
    missing_values=".", 
    filling_values=np.nan  # Replace '.' with NaN
))

k9_median = np.nanmedian(k9_array)
k9_median_plus = k9_median + (k9_median * 1)
print(k9_median_plus)
for i in range(0, k9_array.size):
	if(k9_array[i] < k9_median_plus):
		k9_array[i] = 0
		print("X")
	else:
		k9_array[i] = 1
		print("Y")

#%%

# Step 2: Filter interactions using the binary mask
x = 0
results = []

# Loop through the array
for i in range(array_y_dim):
    chrom_string1, start1, end1 = chrom_coords[i]

    for j in range(array_x_dim):
        chrom_string2, start2, end2 = chrom_coords[j]
        
        # Fetch value and check if it is finite and non-zero
        value = float(chrom_array[i, j])
        if not np.isfinite(value) or value == 0.00:
            continue

        # Compare chromosomes and handle core-chromosome conditions
        if (value >= 3.5) & (chrom_string1 < chrom_string2):
            # Only compute range check once per loop iteration
            # Check if start1 and start2 are within 50kb of the centromere
            if (k9_array[i] == 1 and k9_array[j] == 1):
                results.append([chrom_string1, start1, end1, chrom_string2, start2, end2, "color=color_4_a5", "order1"])
                x += 1
                print(x)
            # Otherwise, add all other interactions
            else:
                results.append([chrom_string1, start1, end1, chrom_string2, start2, end2, "color=color_13_a5", "order0"])
                x += 1



#%%
# Convert to DataFrame
df = pd.DataFrame(results, columns=['Chromosome1', 'Start1', 'End1', 'Chromosome2', 'Start2', 'End2', 'Color', "Order"])

# Sort by order
df = df.sort_values(by="Order")

# Save to file
df.to_csv(f'/Users/iglavinchesk/WP1/WP1.2/config/{name}_5mCcolored_links.txt', sep='\t', index=False, header=False)
# %%
# ALL

# Step 2: Filter interactions using the binary mask
x = 0
results = []

# Loop through the array
for i in range(array_y_dim):
    chrom_string1, start1, end1 = chrom_coords[i]

    for j in range(array_x_dim):
        chrom_string2, start2, end2 = chrom_coords[j]
        
        # Fetch value and check if it is finite and non-zero
        value = float(chrom_array[i, j])
        if not np.isfinite(value) or value == 0.00:
            continue

        # Compare chromosomes and handle core-chromosome conditions
        if (value >= 3.5) & (chrom_string1 < chrom_string2):
            cent_start1, cent_end1 = map(int, cent_array[int(chrom_string1)-1, 1:3])  # Only get start & end
            cent_start2, cent_end2 = map(int, cent_array[int(chrom_string2)-1, 1:3])  # Only get start & end
            # Only compute range check once per loop iteration
            # Check if start1 and start2 are within 50kb of the centromere
            if (cent_start1 - 70000 <= start1 <= cent_start1 + 70000 or cent_end1 - 70000 <= start1 <= cent_end1 + 70000) and \
                (cent_start2 - 70000 <= start2 <= cent_start2 + 70000 or cent_end2 - 70000 <= start2 <= cent_end2 + 70000):
                results.append([chrom_string1, start1, end1, chrom_string2, start2, end2, "color=black_a5", "order1"])
                x += 1
                print(x)
            elif (int(chrom_string1) in core and int(chrom_string2) in accessory) or (int(chrom_string1) in accessory and int(chrom_string2) in core):
                results.append([chrom_string1, start1, end1, chrom_string2, start2, end2, "color=color_5_a5", "order2"])
                x += 1
            elif (k9_array[i] == 1 and k9_array[j] == 1):
                results.append([chrom_string1, start1, end1, chrom_string2, start2, end2, "color=color_4_a5", "order3"])
                x += 1
            # Otherwise, add all other interactions
            else:
                results.append([chrom_string1, start1, end1, chrom_string2, start2, end2, "color=color_13_a5", "order0"])
                x += 1
# %%
# Convert to DataFrame
df = pd.DataFrame(results, columns=['Chromosome1', 'Start1', 'End1', 'Chromosome2', 'Start2', 'End2', 'Color', "Order"])

# Sort by order
df = df.sort_values(by="Order")

# Save to file
df.to_csv(f'/Users/iglavinchesk/WP1/WP1.2/config/{name}_ALLcolored_links.txt', sep='\t', index=False, header=False)

# %%
# ALL-h3k9

# Step 2: Filter interactions using the binary mask
x = 0
results = []

# Loop through the array
for i in range(array_y_dim):
    chrom_string1, start1, end1 = chrom_coords[i]

    for j in range(array_x_dim):
        chrom_string2, start2, end2 = chrom_coords[j]
        
        # Fetch value and check if it is finite and non-zero
        value = float(chrom_array[i, j])
        if not np.isfinite(value) or value == 0.00:
            continue

        # Compare chromosomes and handle core-chromosome conditions
        if (value >= 3.5) & (chrom_string1 < chrom_string2):
            cent_start1, cent_end1 = map(int, cent_array[int(chrom_string1)-1, 1:3])  # Only get start & end
            cent_start2, cent_end2 = map(int, cent_array[int(chrom_string2)-1, 1:3])  # Only get start & end
            # Only compute range check once per loop iteration
            # Check if start1 and start2 are within 50kb of the centromere
            if (cent_start1 - 70000 <= start1 <= cent_start1 + 70000 or cent_end1 - 70000 <= start1 <= cent_end1 + 70000) and \
                (cent_start2 - 70000 <= start2 <= cent_start2 + 70000 or cent_end2 - 70000 <= start2 <= cent_end2 + 70000):
                results.append([chrom_string1, start1, end1, chrom_string2, start2, end2, "color=black_a5", "order1"])
                x += 1
                print(x)
            elif (int(chrom_string1) in core and int(chrom_string2) in accessory) or (int(chrom_string1) in accessory and int(chrom_string2) in core):
                results.append([chrom_string1, start1, end1, chrom_string2, start2, end2, "color=color_5_a5", "order2"])
                x += 1
            #elif (k9_array[i] == 1 and k9_array[j] == 1):
            #    results.append([chrom_string1, start1, end1, chrom_string2, start2, end2, "color=color_4_a5", "order3"])
            #    x += 1
            # Otherwise, add all other interactions
            else:
                results.append([chrom_string1, start1, end1, chrom_string2, start2, end2, "color=color_13_a5", "order0"])
                x += 1
# %%
# Convert to DataFrame
df = pd.DataFrame(results, columns=['Chromosome1', 'Start1', 'End1', 'Chromosome2', 'Start2', 'End2', 'Color', "Order"])

# Sort by order
df = df.sort_values(by="Order")

# Save to file
df.to_csv(f'/Users/iglavinchesk/WP1/WP1.2/config/{name}_cent-acc-colored_links.txt', sep='\t', index=False, header=False)


# %%
