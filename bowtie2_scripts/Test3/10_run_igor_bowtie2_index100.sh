#!/bin/bash -l
# AUTHOR: Yang Chen (yac027@ucsd.edu)
# DATE: 6/30/25
# DESCRIPTION: Igor's Bowtie2 genome-independent read mapping with S. aureus monoculture data against WoLr2 100-partition index

#SBATCH --job-name=bowtie2_igor_index100
#SBATCH --output=/home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_igor_index100_%A_%a.out
#SBATCH --error=/home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_igor_index100_%A_%a.err
#SBATCH --array=0-15
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=96G
#SBATCH --partition=short
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yac027@ucsd.edu

PARTITION_COUNT=100
PARTITION_DIR="/ddn_scratch/y1weng/71_db_splitting/index100"
FASTQ_DIR="/home/yac027/Igor_bowtie2_test/Staph_fastqs"
OUTPUT_DIR="/home/yac027/Igor_bowtie2_test/output_sam_index100_igor"
cd /home/yac027/Igor_bowtie2_test/bowtie2
mkdir -p "$OUTPUT_DIR"

FASTQ_FILES=($FASTQ_DIR/*R1.trimmed.filtered.fastq.gz)
NUM_FILES=${#FASTQ_FILES[@]}
if [[ $SLURM_ARRAY_TASK_ID -ge $NUM_FILES ]]; then
    echo "SLURM task ID exceeds available FASTQ files. Exiting."
    exit 1
fi

R1=${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}
R2="${R1/.R1./.R2.}"
SAMPLE_NAME=$(basename "$R1" .R1.trimmed.filtered.fastq.gz)

echo "Running Igor Bowtie2 on sample $SAMPLE_NAME"

for ((i=1; i<=PARTITION_COUNT; i++)); do
    REF_INDEX="${PARTITION_DIR}/split_${i}"
    OUTFILE="${OUTPUT_DIR}/${SAMPLE_NAME}_index${i}_l100_d.sam"

    ./bowtie2 -p 16 --no-exact-upfront --no-1mm-upfront -x "$REF_INDEX" \
    -1 "$R1" -2 "$R2" --very-sensitive --seed 42 \
    -a --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.05" \
    --no-unal -l 100 -d -S "$OUTFILE"
done

echo "Finished processing $SAMPLE_NAME"

