#!/bin/bash
#SBATCH --job-name=IPO_10kb
#SBATCH --output=hicexplorer_output_%j.log  # %j will be replaced by the job ID
#SBATCH --time=24:00:00                    # Adjust as necessary
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5000
#SBATCH --array=1-4

# Load modules and activate the environment
module load stack
source /cluster/project/gdc/shared/stack/GDCstack.sh
module load python/3.11.6
source /cluster/home/clorrain/miniconda_hic/bin/activate hicexplorer

# Define directories
scripts=/cluster/home/clorrain/Hi-C_analyses/
bamfiles=/cluster/scratch/clorrain/Hi-C/bam/
outdir=/cluster/scratch/clorrain/Hi-C/HicExplorer/H5/
qc_base_dir=/cluster/scratch/clorrain/Hi-C/hicExplorer/H5/
rspath=/cluster/home/clorrain/Hi-C_analyses/

# Define sample-specific restriction cut files and sequences
declare -A restriction_cut_files_ipo323
declare -A restriction_cut_files_ir26b

# Define restriction cut files for IPO323
restriction_cut_files_ipo323["DpnII"]="/cluster/home/clorrain/Hi-C_analyses/IPO323_restrictionCutFile_DpnII.bed"
restriction_cut_files_ipo323["DdeI"]="/cluster/home/clorrain/Hi-C_analyses/IPO323_restrictionCutFile_DdeI.bed"
restriction_cut_files_ipo323["HinfI"]="/cluster/home/clorrain/Hi-C_analyses/IPO323_restrictionCutFile_HinfI.bed"
restriction_cut_files_ipo323["MseI"]="/cluster/home/clorrain/Hi-C_analyses/IPO323_restrictionCutFile_MseI.bed"

restriction_cut_files_ir26b["DpnII"]="/cluster/home/clorrain/Hi-C_analyses/IR26b_restrictionCutFile_DpnII.bed"
restriction_cut_files_ir26b["DdeI"]="/cluster/home/clorrain/Hi-C_analyses/IR26b_restrictionCutFile_DdeI.bed"
restriction_cut_files_ir26b["HinfI"]="/cluster/home/clorrain/Hi-C_analyses/IR26b_restrictionCutFile_HinfI.bed"
restriction_cut_files_ir26b["MseI"]="/cluster/home/clorrain/Hi-C_analyses/IR26b_restrictionCutFile_MseI.bed"

# Build matrices for each sample
#for sample in "${samples[@]}"; do
#    # Extract the base sample name (e.g., IPO323 or IR26b)
#    base_sample=$(echo $sample | cut -d'_' -f1)
#    
#    # Assign the appropriate restriction cut files, sequence, and dangling sequence based on the base sample name
#    if [[ $base_sample == "IPO323" ]]; then
#        restriction_cut_files=("${restriction_cut_files_ipo323[@]}")
#    else
#        restriction_cut_files=("${restriction_cut_files_ir26b[@]}")
#    fi
#    restriction_sequence=restriction_sequences
#    dangling_sequence=restriction_sequences
#    
#    # Define QC directory for the sample
#    qc_dir="${qc_base_dir}${sample}_QC/"
#    
#    # Build the matrix for the sample
#    hicBuildMatrix --samFiles ${bamfiles}${sample}_1.bam ${bamfiles}${sample}_2.bam \
#    -o ${outdir}${sample}_bin1kb.h5 --binSize 1000  --QCfolder ${qc_dir} \
#    --restrictionCutFile "${restriction_cut_files[@]}" --restrictionSequence GATC CT*AG GA*TC TTAA --danglingSequence GATC T*AG A*TC TAA
#done

# Define sample list and select the current sample based on the array task ID
samples=("IPO323_R1" "IPO323_R2" "IPO323+dim2_R1" "IPO323+dim2_R2")
sample=${samples[$SLURM_ARRAY_TASK_ID-1]}

qc_dir=/cluster/scratch/clorrain/Hi-C/hicExplorer/H5/${sample}_QC/

hicBuildMatrix --samFiles ${bamfiles}${sample}_1.bam ${bamfiles}${sample}_2.bam \
-o ${outdir}${sample}_bin10kb.h5 --binSize 10000  --QCfolder ${qc_dir} \
--restrictionCutFile "${restriction_cut_files[@]}" --restrictionSequence GATC CT*AG GA*TC TTAA --danglingSequence GATC T*AG A*TC TAA
