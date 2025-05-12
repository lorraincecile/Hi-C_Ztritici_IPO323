import h5py
import hdf5plugin  # Ensure the plugins are loaded
import numpy as np
from scipy.sparse import csr_matrix
import re

# Function to check the percentage of loci with >= 1000 contacts
def check_contact_threshold(matrix, threshold=1000):
    matrix.setdiag(0)  # Exclude self-interactions by setting the diagonal to zero
    contacts_per_locus = np.array(matrix.sum(axis=1)).flatten()  # Sum contacts per row (locus)
    loci_above_threshold = np.sum(contacts_per_locus >= threshold)
    percentage = (loci_above_threshold / len(contacts_per_locus)) * 100
    return percentage

# Function to process each .h5 file and output results per chromosome
def process_h5_file(file_name):
    try:
        with h5py.File(file_name, 'r') as f:
            # Retrieve chromosome list and other interval data
            chromosomes = f['intervals/chr_list'][:]
            start_list = f['intervals/start_list'][:]
            end_list = f['intervals/end_list'][:]
            
            data = f['matrix/data'][:]
            indices = f['matrix/indices'][:]
            indptr = f['matrix/indptr'][:]
            shape = f['matrix/shape'][:]
            
            # Create the sparse matrix
            contact_matrix = csr_matrix((data, indices, indptr), shape=shape)
    except Exception as e:
        print(f"Error processing {file_name}: {e}")
        return

    # Extract bin size from the filename using regex
    match = re.search(r'_bin(\d+)kb', file_name)
    if match:
        bin_size = int(match.group(1))
    else:
        print(f"Could not extract bin size from {file_name}. Skipping this file.")
        return

    # Store results in a list to write to file later
    results = []

    # Iterate over unique chromosomes
    unique_chromosomes = np.unique(chromosomes)
    for chrom in unique_chromosomes:
        # Identify the rows corresponding to this chromosome
        chrom_indices = np.where(chromosomes == chrom)[0]
        chrom_matrix = contact_matrix[chrom_indices, :][:, chrom_indices]
        
        print(f"Processing {file_name}, Chromosome: {chrom.decode('utf-8')}, Bin size: {bin_size}kb")
        percentage = check_contact_threshold(chrom_matrix)
        
        # Output the result
        results.append(f"{file_name}\t{chrom.decode('utf-8')}\t{bin_size}\t{percentage:.2f}")

    # Output results to a tab-separated file
    output_file = file_name.replace('.h5', '_contact_thresholds.tsv')
    with open(output_file, 'w') as f_out:
        f_out.write("File\tChromosome\tBin_size_kb\tPercentage_loci_above_1000_contacts\n")
        f_out.write("\n".join(results))

# Function to process multiple files
def process_files(files):
    for file in files:
        process_h5_file(file)

# Example usage
if __name__ == "__main__":
    # List of files to process
    files = ["MLP98AG31_H0_bin10kb_final_corrected.h5", "MLP98AG31_H0_bin50kb_final_corrected.h5"]

    # Process files
    process_files(files)