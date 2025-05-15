#!/bin/bash
#SBATCH --job-name=00_RNAseq_trimmomatic
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1g

# --------------------------
# Modules
# --------------------------
module load stack
module load trimmomatic/0.39
module load fastx-toolkit/0.0.14

echo "Modules loaded"
echo "**************"

# --------------------------
# Parameters
# --------------------------
mapping_files=(
    "/cluster/home/iglavinchesk/file_lists/inVitro_RNAseq.txt"
    "/cluster/home/iglavinchesk/file_lists/inPlanta_RNAseq.txt"
)
input_dir="/cluster/scratch/iglavinchesk/05_inPlantaRNAseq/0_raw_data/"
output_dir="/cluster/scratch/iglavinchesk/05_inPlantaRNAseq/01_trimmomatic/"
filtered_dir="/cluster/scratch/iglavinchesk/05_inPlantaRNAseq/02_filtered/"
masked_dir="/cluster/scratch/iglavinchesk/05_inPlantaRNAseq/03_masked/"
adapters="/cluster/home/iglavinchesk/scripts/Tru-Seq3-SE.fa"

trimmomatic_options=(
    "ILLUMINACLIP:${adapters}:2:30:10:2:keepBothReads HEADCROP:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30"
    "ILLUMINACLIP:${adapters}:20:30:15 MINLEN:100"
)

# --------------------------
# Directory Checks
# --------------------------
if [[ ! -d $input_dir ]]; then
    echo "Directory not found: $input_dir"
    exit 1
fi

for dir in "$output_dir" "$filtered_dir" "$masked_dir"; do
    if [[ ! -d $dir ]]; then
        echo "Directory not found. Creating: $dir"
        mkdir -p "$dir"
    fi
done

# --------------------------
# Processing Function
# --------------------------
process_reads() {
    local mapping_file=$1
    local options=$2

    if [[ ! -f $mapping_file ]]; then
        echo "File not found: $mapping_file"
        exit 1
    fi

    while IFS= read -r sample; do
        # Skip lines that start with "#" (header or ignored)
        [[ $sample == \#* ]] && continue

        echo "Processing sample: $sample"

        # Trimming reads
        trimmomatic SE -threads 8 \
            "${input_dir}${sample}.fastq" \
            "${output_dir}${sample}.fastq" \
            $options

        # Quality filtering
        fastq_quality_filter -q 20 -p 80 -v -Q33 \
            -i "${output_dir}${sample}.fastq" \
            -o "${filtered_dir}${sample}-trim-filt.fastq"

        # Masking low-quality bases
        fastq_masker -q 20 -r N -v -Q33 \
            -i "${filtered_dir}${sample}-trim-filt.fastq" \
            -o "${masked_dir}${sample}-trim-filt-mask.fastq"

        echo "Finished processing sample: $sample"
        echo "**************"
    done < "$mapping_file"
}

# --------------------------
# Main Script Execution
# --------------------------
echo "Starting Trimmomatic processing..."

# Process each mapping file with its corresponding options
for i in "${!mapping_files[@]}"; do
    process_reads "${mapping_files[$i]}" "${trimmomatic_options[$i]}"
done

echo "Trimmomatic processing completed."