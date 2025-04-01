#!/bin/bash

# Create output directory if it doesn't exist
# cd output_sam_files_basic/bams
cd output_sam_files_igor_fixed/bams

mkdir -p basic_alignment_stats

# Loop through all .bam files
for bam_file in *.bam; do
    # Extract base filename (without extension)
    base_name=$(basename "$bam_file" .bam)

    # Run samtools flagstat
    samtools flagstat "$bam_file" > "basic_alignment_stats/flagstat_summary_${base_name}.txt"

    # Run samtools stats
    samtools stats "$bam_file" > "basic_alignment_stats/stats_summary_${base_name}.txt"

    echo "Processed: $bam_file"
done

echo "All alignment stats written to basic_alignment_stats/"

