#!/bin/bash
#SBATCH --job-name=02_RNAseq_TEtranscripts
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=1g

# -------------------------
# Directories
# -------------------------
IPO323dir=/cluster/scratch/iglavinchesk/05_inPlantaRNAseq/ref
IPO323_TE_gtf=/cluster/scratch/iglavinchesk/05_inPlantaRNAseq/ref/IPO323_TEs_clean.gtf
IPO323_gene_gtf=/cluster/scratch/iglavinchesk/05_inPlantaRNAseq/ref/IPO323_exons.gtf
STAR_outputDir=/cluster/scratch/iglavinchesk/05_inPlantaRNAseq/02_STAR
outputDir=/cluster/scratch/iglavinchesk/05_inPlantaRNAseq/03_TEtranscripts

# --------------------------
# Variables
# --------------------------
CONDITIONS=("IPO323_A_" "IPO323_B_" "IPO323_C_" "IPO323_D_")
CONTROL="IPO323_WT_invitro_"

# Initialize arrays for condition and control files
declare -A condition_files
control_files=()

# --------------------------
# Assign Files to Variables
# --------------------------
# Populate control files
for replicate in 1 2; do
    control_files+=("${STAR_outputDir}/${CONTROL}${replicate}/${CONTROL}${replicate}.sorted.bam")
done

# Populate condition files for each condition
for condition in "${CONDITIONS[@]}"; do
    for replicate in 1 2; do
        condition_files["$condition"]+="${STAR_outputDir}/${condition}${replicate}/${condition}${replicate}.sorted.bam "
    done
done

# --------------------------
# TEtranscripts Execution
# --------------------------
echo "Starting TEtranscripts analysis..."

for condition in "${CONDITIONS[@]}"; do
    # Define output directory for the condition
    condition_output_dir="${outputDir}/TEtranscripts_${condition}vs_inVitro"
    mkdir -p "$condition_output_dir"

    echo "Running TEtranscripts for condition: $condition"
    echo "Treatment files: ${condition_files[$condition]}"
    echo "Control files: ${control_files[@]}"

    # Run TEtranscripts
    TEtranscripts --treatment ${condition_files[$condition]} \
                  --control ${control_files[@]} \
                  --GTF ${IPO323_gene_gtf} \
                  --TE ${IPO323_TE_gtf} \
                  --format BAM \
                  --mode multi \
                  --outdir "$condition_output_dir"

    echo "Finished TEtranscripts for condition: $condition"
    echo "**************"
done

echo "TEtranscripts analysis completed."