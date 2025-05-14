#!/bin/bash
#SBATCH --job-name=01_ChIPseq_bowtie2_allo
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5g

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
threads=4
bowtie2_options="-k 25 --no-mixed --no-discordant -p $threads"
samtools_options="-o"
allo_options="-seq pe"

# Directories and files
mapping_file_p=/cluster/home/iglavinchesk/file_lists/ChIPseq_paired
mapping_file_s=/cluster/home/iglavinchesk/file_lists/ChIPseq_single
output_dir=/cluster/scratch/iglavinchesk/04_ChIPSeq/02_bowtie2
output_sort=/cluster/scratch/iglavinchesk/04_ChIPSeq/02_bowtie2/sorted
output_allo=/cluster/scratch/iglavinchesk/04_ChIPSeq/03_allo
input_dir=/cluster/scratch/iglavinchesk/04_ChIPSeq/01_trimmomatic/
ref=/cluster/scratch/iglavinchesk/04_ChIPSeq/ref/IPO323

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

# Function to process paired reads
process_paired_reads() {
    echo "Processing paired reads..."
    if [[ ! -f $mapping_file_p ]]; then
        echo "Error: Mapping file for paired reads not found: $mapping_file_p"
        exit 1
    fi

    while IFS= read -r sample; do
        # Skip lines that start with "#" (header or ignored)
        [[ $sample == \#* ]] && continue

        # Check if input files exist
        if [[ ! -f "${input_dir}/${sample}_1_P.fq.gz" || ! -f "${input_dir}/${sample}_2_P.fq.gz" ]]; then
            echo "*ERROR* Input files not found for sample: $sample"
            continue
        fi

        echo "Start ${sample} mapping..."

        # Bowtie2 mapping
        bowtie2 -x ${ref} -1 "${input_dir}/${sample}_1_P.fq.gz" -2 "${input_dir}/${sample}_2_P.fq.gz" \
            -S "${output_dir}/${sample}_bt2.aln.sam" $bowtie2_options

        echo "Sorting..."
        samtools collate $samtools_options "${output_sort}/${sample}_bt2.sort.sam" "${output_dir}/${sample}_bt2.aln.sam"

        echo "Remapping MMRs..."
        allo "${output_sort}/${sample}_bt2.sort.sam" $allo_options -o "${output_allo}/${sample}_allo"

        echo "Finish processing ${sample}"
        echo "**************"

    done < "$mapping_file_p"
}

# Function to process single reads
process_single_reads() {
    echo "Processing single reads..."
    if [[ ! -f $mapping_file_s ]]; then
        echo "Error: Mapping file for single reads not found: $mapping_file_s"
        exit 1
    fi

    while IFS= read -r sample; do
        # Skip lines that start with "#" (header or ignored)
        [[ $sample == \#* ]] && continue

        # Check if input file exists
        if [[ ! -f "${input_dir}/${sample}.fq.gz" ]]; then
            echo "*ERROR* Input file not found for sample: $sample"
            continue
        fi

        echo "Start ${sample} mapping..."

        # Bowtie2 mapping
        bowtie2 -x ${ref} -q "${input_dir}/${sample}.fq.gz" \
            -S "${output_dir}/${sample}_bt2.aln.sam" $bowtie2_options

        echo "Sorting..."
        samtools collate $samtools_options "${output_sort}/${sample}_bt2.sort.sam" "${output_dir}/${sample}_bt2.aln.sam"

        echo "Remapping MMRs..."
        allo "${output_sort}/${sample}_bt2.sort.sam" $allo_options -o "${output_allo}/${sample}_allo"

        echo "Finish processing ${sample}"
        echo "**************"

    done < "$mapping_file_s"
}

# --------------------------
# Main Script Execution
# --------------------------
echo "Starting Bowtie2 and Allo processing..."

# Check and create necessary directories
check_and_create_dir "$output_dir"
check_and_create_dir "$output_sort"
check_and_create_dir "$output_allo"

# Process paired and single reads
process_paired_reads
process_single_reads

echo "Bowtie2 and Allo processing completed."