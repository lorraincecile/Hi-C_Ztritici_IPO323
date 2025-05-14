#!/bin/bash
#SBATCH --job-name=02_ChIPseq_sort_index_bam
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=750m

# --------------------------
# Modules
# --------------------------
module load stack
module load samtools
source /cluster/project/gdc/shared/stack/GDCstack.sh
module load stack/2024-06 gcc/12.2.0
module load python
module load bowtie2/2.5.1-u2j3omo
source ${HOME}/packages/allo/bin/activate

echo "Modules loaded"
echo "**************"

# --------------------------
# Parameters
# --------------------------
bam_quality=30
samtools_options="-bS -q $bam_quality"
mapping_file_p=/cluster/home/iglavinchesk/file_lists/ChIPseq_paired
mapping_file_s=/cluster/home/iglavinchesk/file_lists/ChIPseq_single
output_allo=/cluster/scratch/iglavinchesk/04_ChIPSeq/03_allo
sorted_dir="${output_allo}/sorted"

# --------------------------
# Functions
# --------------------------

# Function to check and create directories
check_and_create_dir() {
    if [[ ! -d $1 ]]; then
        echo "Directory not found. Creating: $1"
        mkdir -p "$1"
    fi
}

# Function to process BAM files
process_bam_files() {
    local mapping_file=$1
    local seq_type=$2

    echo "Processing $seq_type reads..."
    if [[ ! -f $mapping_file ]]; then
        echo "Error: Mapping file not found: $mapping_file"
        exit 1
    fi

    while IFS= read -r sample; do
        # Skip lines that start with "#" (header or ignored)
        [[ $sample == \#* ]] && continue

        # Check if input file exists
        if [[ ! -f "${output_allo}/${sample}_allo" ]]; then
            echo "*ERROR* Input file not found for sample: $sample"
            continue
        fi

        echo "Processing sample: $sample ($seq_type reads)"

        # Convert SAM to BAM
        echo "Converting SAM to BAM..."
        samtools view $samtools_options "${output_allo}/${sample}_allo" > "${output_allo}/${sample}_allo.q${bam_quality}.bam"

        # Remove original SAM file
        rm -f "${output_allo}/${sample}_allo"

        # Sort BAM file
        echo "Sorting BAM..."
        samtools sort "${output_allo}/${sample}_allo.q${bam_quality}.bam" > "${sorted_dir}/${sample}_sorted.allo.q${bam_quality}.bam"

        # Index BAM file
        echo "Indexing BAM..."
        samtools index "${sorted_dir}/${sample}_sorted.allo.q${bam_quality}.bam"

        echo "Finished processing sample: $sample"
        echo "**************"

    done < "$mapping_file"
}

# --------------------------
# Main Script Execution
# --------------------------
echo "Starting BAM sorting and indexing..."

# Check and create necessary directories
check_and_create_dir "$output_allo"
check_and_create_dir "$sorted_dir"

# Process paired and single reads
process_bam_files "$mapping_file_p" "paired"
process_bam_files "$mapping_file_s" "single"

echo "BAM sorting and indexing completed."