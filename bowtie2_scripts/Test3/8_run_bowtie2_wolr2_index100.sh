#!/bin/bash -l
# AUTHOR: Yang Chen (yac027@ucsd.edu)
# DATE: 6/30/25
# DESCRIPTION: Bowtie2 alignment of S. aureus monoculture data against WoLr2 100-partition index

#SBATCH --job-name=bowtie2_wolr2_index100
#SBATCH --output=/home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_wolr2_index100_%A_%a.out
#SBATCH --error=/home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_wolr2_index100_%A_%a.err
#SBATCH --array=0-15
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yac027@ucsd.edu

# Activate Bowtie2 environment
source activate bowtie2

PARTITION_COUNT=100
PARTITION_DIR="/ddn_scratch/y1weng/71_db_splitting/index100"
FASTQ_DIR="/home/yac027/Igor_bowtie2_test/Staph_fastqs"
OUTPUT_DIR="/home/yac027/Igor_bowtie2_test/output_sam_index100"
LOG_DIR="/home/yac027/Igor_bowtie2_test/slurm_out"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

which bowtie2 | tee -a "$LOG_DIR/bowtie2_wolr2_index100_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

FASTQ_FILES=($FASTQ_DIR/*R1.trimmed.filtered.fastq.gz)
NUM_FILES=${#FASTQ_FILES[@]}

if [[ $SLURM_ARRAY_TASK_ID -ge $NUM_FILES ]]; then
    echo "SLURM task ID exceeds input file count." | tee -a "$LOG_DIR/bowtie2_wolr2_index100_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
    exit 1
fi

FASTQ1=${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}
FASTQ2="${FASTQ1/.R1./.R2.}"
SAMPLE_NAME=$(basename "$FASTQ1" .R1.trimmed.filtered.fastq.gz)

echo "Processing sample: $SAMPLE_NAME" | tee -a "$LOG_DIR/bowtie2_wolr2_index100_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

for ((i=1; i<=PARTITION_COUNT; i++)); do
    PART_INDEX="${PARTITION_DIR}/split_${i}"
    OUTFILE="${OUTPUT_DIR}/${SAMPLE_NAME}_index${i}.sam"

    echo "Mapping to index partition $i..." | tee -a "$LOG_DIR/bowtie2_wolr2_index100_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

    bowtie2 -p 16 \
      --no-exact-upfront --no-1mm-upfront \
      --very-sensitive \
      --seed 42 \
      -k 16 \
      --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" \
      --score-min "L,0,-0.05" \
      --no-unal \
      -x "$PART_INDEX" \
      -1 "$FASTQ1" -2 "$FASTQ2" \
      -S "$OUTFILE" 2>> "$LOG_DIR/bowtie2_wolr2_index100_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

    echo "Finished index $i" | tee -a "$LOG_DIR/bowtie2_wolr2_index100_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
done

echo "Completed all partitions for $SAMPLE_NAME" | tee -a "$LOG_DIR/bowtie2_wolr2_index100_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

