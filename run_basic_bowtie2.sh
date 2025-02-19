#!/bin/bash -l
# AUTHOR: Yang Chen (yac027@ucsd.edu)
# DATE: 2/10/25
# DESCRIPTION: Minimal Bowtie2 test with S. aureus monoculture data

#SBATCH --job-name=bowtie2_basic
#SBATCH --output=/home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_basic_%A_%a.out
#SBATCH --error=/home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_basic_%A_%a.err
#SBATCH --array=0-15
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=48G
#SBATCH --partition=short

# Activate Bowtie2 environment
source activate bowtie2

# Set relevant paths
REF_INDEX="/home/yac027/Igor_bowtie2_test/S.aureus_reference/index_files/GCF_001456215.1_bowtie2_index"
FASTQ_DIR="/home/yac027/Igor_bowtie2_test/Staph_fastqs"
OUTPUT_DIR="/home/yac027/Igor_bowtie2_test/output_sam_files_basic"
mkdir -p "$OUTPUT_DIR"

# Log Bowtie2 path
which bowtie2 | tee -a /home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_basic_%A_%a.log

# Find FASTQ files
FASTQ_FILES=($FASTQ_DIR/*R1.trimmed.filtered.fastq.gz)
NUM_FILES=${#FASTQ_FILES[@]}

# Check if task ID is within range
if [[ $SLURM_ARRAY_TASK_ID -ge $NUM_FILES ]]; then
    echo "SLURM task ID exceeds available FASTQ files. Exiting." | tee -a /home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_basic_%A_%a.log
    exit 1
fi

# Assign FASTQ files
FASTQ1=${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}
FASTQ2="${FASTQ1/.R1./.R2.}"
SAMPLE_NAME=$(basename "$FASTQ1" .R1.trimmed.filtered.fastq.gz)

# Debugging check
echo "Processing sample: $SAMPLE_NAME" | tee -a /home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_basic_%A_%a.log
echo "FASTQ1: $FASTQ1" | tee -a /home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_basic_%A_%a.log
echo "FASTQ2: $FASTQ2" | tee -a /home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_basic_%A_%a.log

# Run Bowtie2 with minimal parameters
bowtie2 -x "$REF_INDEX" \
  -1 "$FASTQ1" \
  -2 "$FASTQ2" \
  -S "$OUTPUT_DIR/${SAMPLE_NAME}.sam" \
2>&1 | tee -a /home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_basic_%A_%a.log

echo "Finished processing $SAMPLE_NAME" | tee -a /home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_basic_%A_%a.log

