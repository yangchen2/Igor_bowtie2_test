#!/bin/bash -l
# AUTHOR: Yang Chen (yac027@ucsd.edu)
# DATE: 4/17/25
# DESCRIPTION: Test Igor's branch of Bowtie2 genome-independent read mapping with S. aureus monoculture data against WoLr2 database

#SBATCH --job-name=bowtie2_igor_wolr2
#SBATCH --output=/home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_igor_wolr2_%A_%a.out
#SBATCH --error=/home/yac027/Igor_bowtie2_test/slurm_out/bowtie2_igor_wolr2_%A_%a.err
#SBATCH --array=0-15
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=96G
#SBATCH --partition=short

echo "Setting relevant paths"
REF_INDEX="/scratch/qp-woltka/WoLr2/WoLr2"
FASTQ_DIR="/home/yac027/Igor_bowtie2_test/Staph_fastqs"
OUTPUT_DIR="/home/yac027/Igor_bowtie2_test/output_sam_files_igor_wolr2"
mkdir -p $OUTPUT_DIR

echo "Finding FASTQ files"
FASTQ_FILES=($FASTQ_DIR/*R1.trimmed.filtered.fastq.gz)
NUM_FILES=${#FASTQ_FILES[@]}

if [[ $SLURM_ARRAY_TASK_ID -ge $NUM_FILES ]]; then
    echo "SLURM task ID exceeds available FASTQ files. Exiting."
    exit 1
fi

R1=${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}
R2="${R1/.R1./.R2.}"

# Debugging check
echo "R1: $R1"
echo "R2: $R2"

SAMPLE_NAME=$(basename "$R1" .R1.trimmed.filtered.fastq.gz)

echo "Running Bowtie2 on sample $SAMPLE_NAME"
cd /home/yac027/Igor_bowtie2_test/bowtie2
./bowtie2 -p 16 --no-exact-upfront --no-1mm-upfront -x $REF_INDEX \
-1 "$R1" -2 "$R2" --very-sensitive --seed 42 \
-a --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.05" \
--no-unal -l 100 -d -S "$OUTPUT_DIR/${SAMPLE_NAME}_l100_d.sam"

echo "Finished processing $SAMPLE_NAME"

