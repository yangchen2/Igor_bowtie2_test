#!/bin/bash -l
# AUTHOR: Yang Chen (yac027@ucsd.edu)
# DATE: 4/17/25
# DESCRIPTION: Bowtie2 alignment with -a flag to retain all alignments against WoLr2 database

#SBATCH --job-name=bowtie2_wolr2_all
#SBATCH --output=/home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_wolr2_all_%A_%a.out
#SBATCH --error=/home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_wolr2_all_%A_%a.err
#SBATCH --array=0-15
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=96G
#SBATCH --partition=short

# Activate Bowtie2 environment
source activate bowtie2

# Set paths
REF_INDEX="/scratch/qp-woltka/WoLr2/WoLr2"
FASTQ_DIR="/home/yac027/Igor_bowtie2_test/Staph_fastqs"
OUTPUT_DIR="/home/yac027/Igor_bowtie2_test/output_sam_files_wolr2_all"
mkdir -p "$OUTPUT_DIR"

# Log Bowtie2 path
which bowtie2 | tee -a /home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_wolr2_all_%A_%a.log

# Get list of FASTQ R1 files
FASTQ_FILES=($FASTQ_DIR/*R1.trimmed.filtered.fastq.gz)
NUM_FILES=${#FASTQ_FILES[@]}

# Ensure valid task ID
if [[ $SLURM_ARRAY_TASK_ID -ge $NUM_FILES ]]; then
    echo "SLURM task ID exceeds available FASTQ files. Exiting." | tee -a /home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_wolr2_all_%A_%a.log
    exit 1
fi

# Get input files and sample name
FASTQ1=${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}
FASTQ2="${FASTQ1/.R1./.R2.}"
SAMPLE_NAME=$(basename "$FASTQ1" .R1.trimmed.filtered.fastq.gz)

echo "Processing sample: $SAMPLE_NAME" | tee -a /home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_wolr2_all_%A_%a.log
echo "FASTQ1: $FASTQ1" | tee -a /home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_wolr2_all_%A_%a.log
echo "FASTQ2: $FASTQ2" | tee -a /home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_wolr2_all_%A_%a.log

# Run Bowtie2 alignment with -a flag against WoLr2 index
bowtie2 -p 16 \
  --no-exact-upfront --no-1mm-upfront \
  --very-sensitive \
  --seed 42 \
  -a \
  --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" \
  --score-min "L,0,-0.05" \
  --no-unal \
  -x "$REF_INDEX" \
  -1 "$FASTQ1" -2 "$FASTQ2" \
  -S "$OUTPUT_DIR/${SAMPLE_NAME}.sam"

echo "Finished processing $SAMPLE_NAME" | tee -a /home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_wolr2_all_%A_%a.log

