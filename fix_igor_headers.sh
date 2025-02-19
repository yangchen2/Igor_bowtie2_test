#!/bin/bash -l
#SBATCH --job-name=fix_igor_sam_headers
#SBATCH --output=slurm_out/fix_igor_sam_headers_%j.out
#SBATCH --error=slurm_out/fix_igor_sam_headers_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=short

# Define directories
BASIC_DIR="output_sam_files_basic"
IGOR_DIR="output_sam_files_igor"
FIXED_DIR="output_sam_files_igor_fixed"

# Create an output directory for fixed SAM files
mkdir -p ${FIXED_DIR}

# Loop through each SAM file in the Igor directory
for igor_sam in ${IGOR_DIR}/*.sam; do
    # Extract the sample prefix (remove directory path and _l100_d suffix)
    sample_name=$(basename "$igor_sam" | sed 's/_l100_d.sam//')

    # Define the corresponding basic SAM file
    basic_sam="${BASIC_DIR}/${sample_name}.sam"

    # Check if the basic SAM file exists
    if [[ -f "$basic_sam" ]]; then
        # Extract the header from the basic SAM file
        samtools view -H "$basic_sam" > "${FIXED_DIR}/${sample_name}_fixed.sam"

        # Append the reads from the Igor SAM file
        cat "$igor_sam" >> "${FIXED_DIR}/${sample_name}_fixed.sam"

        echo "Fixed: ${sample_name}_fixed.sam"
    else
        echo "WARNING: No matching basic SAM file found for ${igor_sam}"
    fi
done

echo "Header fixing completed!"

