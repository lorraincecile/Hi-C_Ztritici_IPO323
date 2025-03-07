#!/bin/bash
#SBATCH --mem-per-cpu=4000 
#SBATCH --cpus-per-task=3 
#SBATCH --time=4:00:00 


source $1 
#IDX=${LSB_JOBINDEX}
IDX=${SLURM_ARRAY_TASK_ID}

# Load modules
# ------------
#module load gdc java/1.8.0_31 gatk/4.1.2.0
source /cluster/project/gdc/shared/stack/GDCstack.sh
module load trimmomatic/0.39
module load bwa
module load samtools 

#Directories
#refgenome=/cluster/home/clorrain/ref_genomes/IPO323_nomt_no18.fa
refgenome=/cluster/home/clorrain/ref_genomes/IR01_26b.fa
read_dir=/cluster/scratch/clorrain/Hi-C/IR26b-dim2_R1/
outdir=/cluster/scratch/clorrain/Hi-C/IR26b-dim2_R1/
sample=IR26b-dim2_R1
adapters=adapters=/cluster/home/clorrain/adapters_novogene.fa

#Trimming
# --------
trimmomatic PE \
  ${read_dir}${sample}_1.fq.gz ${read_dir}${sample}_2.fq.gz \
  ${outdir}${sample}_1P.fq.gz ${outdir}${sample}_1U.fq.gz \
  ${outdir}${sample}_2P.fq.gz ${outdir}${sample}_2U.fq.gz \
  ILLUMINACLIP:${adapters}:2:30:10 \
  LEADING:15 TRAILING:15 SLIDINGWINDOW:5:15 MINLEN:50


#Mapping
# --------

#read1
bwa mem -A1 -B4 -E50 -L0 ${refgenome}\
 ${read_dir}${sample}_1P.fq.gz\
 2>>${outdir}${sample}_1.log > ${outdir}${sample}_1.sam
samtools view -Shb ${outdir}${sample}_1.sam  > ${outdir}${sample}_notmt_1.bam

#read2
bwa mem -A1 -B4  -E50 -L0 ${refgenome}\
 ${read_dir}${sample}_2P.fq.gz 2>>${outdir}${sample}_2.log > ${outdir}${sample}_2.sam
samtools view -Shb ${outdir}${sample}_2.sam  > ${outdir}${sample}_2.bam


