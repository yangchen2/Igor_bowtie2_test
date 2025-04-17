#!/bin/bash -l
#SBATCH --job-name=sam_to_bam
#SBATCH --output=slurm_out/sam_to_bam_%j.out
#SBATCH --error=slurm_out/sam_to_bam_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=short

# Load necessary modules (if required, adjust based on your cluster)
source activate samtools

# Define directories
BASIC_DIR="output_sam_files_basic"
# IGOR_DIR="output_sam_files_igor"
IGOR_DIR="output_sam_files_igor_fixed"

# Create an output directory for BAM files
mkdir -p ${BASIC_DIR}/bams
mkdir -p ${IGOR_DIR}/bams

# Function to convert SAM to sorted and indexed BAM
convert_sam_to_bam() {
    sam_file=$1
    output_dir=$2
    base_name=$(basename "$sam_file" .sam)
    
    # Convert to BAM
    samtools view -b -o "$output_dir/${base_name}.bam" "$sam_file"
    
    # Sort BAM
    samtools sort -o "$output_dir/${base_name}_sorted.bam" "$output_dir/${base_name}.bam"
    
    # Index BAM
    samtools index "$output_dir/${base_name}_sorted.bam"
    
    # Remove intermediate BAM
    rm "$output_dir/${base_name}.bam"
}

export -f convert_sam_to_bam

# Convert all SAM files in the basic directory
find $BASIC_DIR -name "*.sam" | parallel -j 4 convert_sam_to_bam {} ${BASIC_DIR}/bams

# Convert all SAM files in the Igor directory
find $IGOR_DIR -name "*.sam" | parallel -j 4 convert_sam_to_bam {} ${IGOR_DIR}/bams

echo "SAM to BAM conversion completed!"

