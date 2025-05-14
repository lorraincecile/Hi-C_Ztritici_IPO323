#!/bin/bash
#SBATCH --job-name=04_normalise2bedGraph
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=200m

# --------------------------
# Modules
# --------------------------
module load stack
module load bedtools2/2.31.0

echo "Modules loaded"
echo "**************"

# --------------------------
# Parameters
# --------------------------
mapping_file=/cluster/home/iglavinchesk/file_lists/ChIPseq_all
input_dir=/cluster/scratch/iglavinchesk/04_ChIPSeq/04_markdup
final_dir=/cluster/scratch/iglavinchesk/04_ChIPSeq/5_normalisation1x/final/10bp
bigwig_dir=/cluster/scratch/iglavinchesk/04_ChIPSeq/5_normalisation1x/bigwig
output_dir=/cluster/scratch/iglavinchesk/04_ChIPSeq/06_macs3/allo.q30.1x
effectiveGenomeSize=39111209
bamCoverage_options="--normalizeUsing RPGC --effectiveGenomeSize $effectiveGenomeSize --binSize 10 --extendReads 150 --smoothLength 150 --centerReads --outFileFormat bedgraph"

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
        if [[ ! -f "${input_dir}/${sample}_markdup.allo.q30.bam" ]]; then
            echo "*ERROR* Input file not found for sample: $sample"
            continue
        fi

        echo "Creating normalized BigWig for ${sample}..."

        # Generate normalized BigWig
        bamCoverage --bam "${input_dir}/${sample}_markdup.allo.q30.bam" \
            -o "${bigwig_dir}/${sample}.allo.q30.1x.bedGraph" \
            -p 8 $bamCoverage_options

        # Sort the bedGraph
        sort -k1,1 -k2,2n "${bigwig_dir}/${sample}.allo.q30.1x.bedGraph" > "${bigwig_dir}/${sample}.allo.q30.1x.sorted.bedGraph"

        echo "Finished processing ${sample}"
        echo "**************"

    done < "$mapping_file"
}

# Function to calculate mean coverage
calculate_mean() {
    local histone=$1
    local output_file=$2
    shift 2
    local input_files=("$@")

    echo "Calculating mean for $histone..."
    bedtools unionbedg -i "${input_files[@]}" | \
    awk '{sum=0; count=0; for(i=4; i<=NF; i++) { if($i != ".") { sum+=$i; count+=1; } } if(count>0) {printf "%s\t%s\t%s\t%.6f\n", $1, $2, $3, sum/count; } else { printf "%s\t%s\t%s\t.\n", $1, $2, $3; }}' > "$output_file"
    echo "Mean calculated for $histone: $output_file"
}

# --------------------------
# Main Script Execution
# --------------------------
echo "Starting normalization and mean calculation..."

# Check and create necessary directories
check_and_create_dir "$final_dir"
check_and_create_dir "$bigwig_dir"
check_and_create_dir "$output_dir"

# Process samples
process_samples

# Calculate means for histone marks
calculate_mean "H3K36me3" "${final_dir}/H3K36me3_merged_mean.RPGC.bedgraph" \
    "${bigwig_dir}/H3K36me3_rep1_2.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H3K36me3_rep2_2.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H3K36me3_rep3.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H3K36me3_rep4.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H3K36me3_rep5.allo.q30.1x.sorted.bedGraph"

calculate_mean "CenH3" "${final_dir}/CenH3_merged_mean.RPGC.bedgraph" \
    "${bigwig_dir}/CenH3_rep1.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/CenH3_rep2.allo.q30.1x.sorted.bedGraph"

calculate_mean "H4K20me3" "${final_dir}/H4K20me3_merged_mean.RPGC.bedgraph" \
    "${bigwig_dir}/H4K20me3_rep1_2.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H4K20me3_rep2_2.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H4K20me3_rep3.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H4K20me3_rep4.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H4K20me3_rep5.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H4K20me3_rep6.allo.q30.1x.sorted.bedGraph"

calculate_mean "H4K20me1" "${final_dir}/H4K20me1_merged_mean.RPGC.bedgraph" \
    "${bigwig_dir}/H4K20me1_rep1.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H4K20me1_rep2.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H4K20me1_rep3.allo.q30.1x.sorted.bedGraph"

calculate_mean "H3K27me3" "${final_dir}/H3K27me3_merged_mean.RPGC.bedgraph" \
    "${bigwig_dir}/H3K27me3_rep1.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H3K27me3_rep2.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H3K27me3_rep3.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H3K27me3_rep4.allo.q30.1x.sorted.bedGraph"

calculate_mean "H3K27me2" "${final_dir}/H3K27me2_merged_mean.RPGC.bedgraph" \
    "${bigwig_dir}/H3K27me2_rep1.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H3K27me2_rep2.allo.q30.1x.sorted.bedGraph"

calculate_mean "H3K9me3" "${final_dir}/H3K9me3_merged_mean.RPGC.bedgraph" \
    "${bigwig_dir}/H3K9me3_rep1.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H3K9me3_rep2.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H3K9me3_rep3.allo.q30.1x.sorted.bedGraph"

calculate_mean "H3K4me2" "${final_dir}/H3K4me2_merged_mean.RPGC.bedgraph" \
    "${bigwig_dir}/H3K4me2_rep1.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H3K4me2_rep2.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H3K4me2_rep3.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H3K4me2_rep4.allo.q30.1x.sorted.bedGraph"

calculate_mean "H3K9me2" "${final_dir}/H3K9me2_merged_mean.RPGC.bedgraph" \
    "${bigwig_dir}/H3K9me2_rep1.allo.q30.1x.sorted.bedGraph" \
    "${bigwig_dir}/H3K9me2_rep2.allo.q30.1x.sorted.bedGraph"

echo "Normalization and mean calculation completed."