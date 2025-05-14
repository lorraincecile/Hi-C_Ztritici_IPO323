#!/bin/bash
#SBATCH --job-name=00_ChIPseq_trimmomatic
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2g

# --------------------------
# Modules
# --------------------------
module load stack
module load trimmomatic/0.39

echo "00 Modules loaded"

# --------------------------
# Parameters
# --------------------------
threads=4
trimmomatic_options="ILLUMINACLIP:TruSeq3-SE.fa:2:30:10:2:keepBothReads HEADCROP:5 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30"
single_trimmomatic_options="ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30"
dry_run=false  # Set to true for testing without running trimmomatic

# Directories and files
mapping_file_p=/cluster/home/iglavinchesk/file_lists/ChIPseq_paired
mapping_file_s=/cluster/home/iglavinchesk/file_lists/ChIPseq_single
raw_data_dir=/cluster/scratch/iglavinchesk/04_ChIPSeq/01_trimmomatic/00_raw_data
input_dir=/cluster/scratch/iglavinchesk/04_ChIPSeq/00_raw_data/
output_dir=/cluster/scratch/iglavinchesk/04_ChIPSeq/01_trimmomatic/
adapters=/cluster/work/gdc/shared/p958/0_raw_data/2_RNAseq/1_scripts/adapters_novogene.fa

# --------------------------
# Check Directories
# --------------------------
if [[ ! -d $input_dir ]]; then
    echo "Error: Input directory not found: $input_dir"
    exit 1
fi

if [[ ! -d $output_dir ]]; then
    echo "Output directory not found. Creating: $output_dir"
    mkdir -p "$output_dir"
fi

# --------------------------
# Function to Process Paired Reads
# --------------------------
process_paired_reads() {
    echo "Processing paired reads..."
    if [[ ! -f $mapping_file_p ]]; then
        echo "Error: Mapping file for paired reads not found: $mapping_file_p"
        exit 1
    fi

    while IFS= read -r sample; do
        # Skip lines that start with "#" (header or ignored)
        [[ $sample == \#* ]] && continue

        echo "Processing sample: $sample (paired reads)"
        if [[ $dry_run == true ]]; then
            echo "Dry run: Skipping trimmomatic for $sample"
            continue
        fi

        trimmomatic PE -threads $threads \
            "${input_dir}${sample}_R1.fastq" \
            "${input_dir}${sample}_R2.fastq" \
            "${output_dir}${sample}_1_P.fq.gz" \
            "${output_dir}${sample}_1_U.fq.gz" \
            "${output_dir}${sample}_2_P.fq.gz" \
            "${output_dir}${sample}_2_U.fq.gz" \
            $trimmomatic_options

        # Remove unpaired reads
        rm -f "${output_dir}${sample}_1_U.fq.gz" "${output_dir}${sample}_2_U.fq.gz"
    done < "$mapping_file_p"
}

# --------------------------
# Function to Process Single Reads
# --------------------------
process_single_reads() {
    echo "Processing single reads..."
    if [[ ! -f $mapping_file_s ]]; then
        echo "Error: Mapping file for single reads not found: $mapping_file_s"
        exit 1
    fi

    while IFS= read -r sample; do
        # Skip lines that start with "#" (header or ignored)
        [[ $sample == \#* ]] && continue

        echo "Processing sample: $sample (single reads)"
        if [[ $dry_run == true ]]; then
            echo "Dry run: Skipping trimmomatic for $sample"
            continue
        fi

        trimmomatic SE -threads $threads \
            "${input_dir}${sample}.fastq" \
            "${output_dir}${sample}.fq.gz" \
            $single_trimmomatic_options
    done < "$mapping_file_s"
}

# --------------------------
# Main Script Execution
# --------------------------
echo "Starting Trimmomatic processing..."
process_paired_reads
process_single_reads
echo "Trimmomatic processing completed."