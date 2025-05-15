#!/bin/bash
#SBATCH --job-name=01_RNAseq_STARalignment
#SBATCH --time=9:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=1g

# -------------------------
# Directories
# -------------------------
IPO323dir=/cluster/scratch/iglavinchesk/05_inPlantaRNAseq/ref
IPO323_TE_gtf=/cluster/scratch/iglavinchesk/05_inPlantaRNAseq/ref/IPO323_TEs_clean.gtf
IPO323_gene_gtf=/cluster/scratch/iglavinchesk/05_inPlantaRNAseq/ref/IPO323_exons.gtf
inputDir=/cluster/scratch/iglavinchesk/05_inPlantaRNAseq/01_trimmomatic
outputDir=/cluster/scratch/iglavinchesk/05_inPlantaRNAseq/02_STAR
IPO323list=/cluster/home/iglavinchesk/file_lists/RNAseq.txt

# --------------------------
# Modules
# --------------------------
module load stack
module load star/2.7.10b
module load samtools/1.16.1

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

# Function to align reads with STAR
align_reads() {
    local sample=$1

    echo "Aligning ${sample}..."

    # Create output directory for the sample
    check_and_create_dir "${outputDir}/${sample}"

    # Run STAR alignment
    STAR \
        --genomeDir "${IPO323dir}" \
        --readFilesIn "${inputDir}/${sample}.fastq" \
        --outFileNamePrefix "${outputDir}/${sample}/${sample}." \
        --outFilterMultimapNmax 200 \
        --winAnchorMultimapNmax 200 \
        --outSAMtype BAM Unsorted \
        --outReadsUnmapped None \
        --runThreadN 30 \
        --limitBAMsortRAM 40000000 \
        --genomeSAindexNbases 14

    echo "Finished aligning ${sample}"

    # Sort the resulting BAM file
    echo "Sorting BAM file for ${sample}..."
    samtools sort -@ 30 -o "${outputDir}/${sample}/${sample}.sorted.bam" "${outputDir}/${sample}/${sample}.Aligned.out.bam"

    # Index the sorted BAM file
    echo "Indexing BAM file for ${sample}..."
    samtools index "${outputDir}/${sample}/${sample}.sorted.bam"

    # Remove the unsorted BAM file to save space
    rm -f "${outputDir}/${sample}/${sample}.Aligned.out.bam"

    echo "Finished processing ${sample}"
    echo "**************"
}

# --------------------------
# Main Script Execution
# --------------------------
echo "Starting STAR alignment..."

# Check and create necessary directories
check_and_create_dir "$inputDir"
check_and_create_dir "$outputDir"

# Process each sample in the list
while IFS= read -r sample; do
    # Skip lines that start with "#" (header or ignored)
    [[ $sample == \#* ]] && continue

    # Align reads for the sample
    align_reads "$sample"
done < "$IPO323list"

echo "STAR alignment completed."