#!/bin/bash
#SBATCH --job-name=03_ChIPseq_picard_markdup
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=20g

# --------------------------
# Modules
# --------------------------
module load stack
module load picard/2.26.2
module load samtools
echo "Modules loaded"
echo "**************"

# --------------------------
# Parameters
# --------------------------
mapping_file=/cluster/scratch/iglavinchesk/04_ChIPSeq/extra.txt
input_dir=/cluster/scratch/iglavinchesk/04_ChIPSeq/07_allo/sorted
output_dir=/cluster/scratch/iglavinchesk/04_ChIPSeq/08_markdup
metrics_dir=/cluster/scratch/iglavinchesk/04_ChIPSeq/08_markdup/metrics
picard_options="VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false"

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

# Function to process samples
process_samples() {
    echo "Processing samples..."
    if [[ ! -f $mapping_file ]]; then
        echo "Error: Mapping file not found: $mapping_file"
        exit 1
    fi

    while IFS= read -r sample; do
        # Skip lines that start with "#" (header or ignored)
        [[ $sample == \#* ]] && continue

        # Check if input file exists
        if [[ ! -f "${input_dir}/${sample}_sorted.allo.q30.bam" ]]; then
            echo "*ERROR* Input file not found for sample: $sample"
            continue
        fi

        echo "Start processing ${sample}..."

        # Mark duplicates using Picard
        picard MarkDuplicates \
            I="${input_dir}/${sample}_sorted.allo.q30.bam" \
            O="${output_dir}/${sample}_markdup.allo.q30.bam" \
            M="${metrics_dir}/${sample}_markdup.allo.q30.txt" \
            $picard_options

        # Index the BAM file
        echo "Indexing BAM..."
        samtools index "${output_dir}/${sample}_markdup.allo.q30.bam"

        echo "Finished processing ${sample}"
        echo "**************"

    done < "$mapping_file"
}

# --------------------------
# Main Script Execution
# --------------------------
echo "Starting Picard MarkDuplicates processing..."

# Check and create necessary directories
check_and_create_dir "$output_dir"
check_and_create_dir "$metrics_dir"

# Process samples
process_samples

echo "Picard MarkDuplicates processing completed."