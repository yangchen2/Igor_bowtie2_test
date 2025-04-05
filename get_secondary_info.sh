#!/bin/bash

# Input SAM files
SAM_BASIC="output_sam_files_basic/NC.Staphy2.1A.sam"
SAM_ALL="output_sam_files_basic_all-reads/NC.Staphy2.1A.sam"
SAM_IGOR=$(ls output_sam_files_igor/NC.Staphy2.1A*.sam)  # wildcard handles suffixes

# Output summary files
OUT_BASIC="secondary_summary_bt2-basic.txt"
OUT_ALL="secondary_summary_bt2-basic_all-reads.txt"
OUT_IGOR="secondary_summary_bt2-igor.txt"

# Function to extract secondary alignments summary
extract_secondary_info() {
  local infile=$1
  local outfile=$2

  samtools view -f 256 "$infile" | awk -F'\t' '
  {
    read=$1; ref=$3; pos=$4; score="NA";
    for (i=12; i<=NF; i++) {
      if ($i ~ /^AS:i:/) {
        split($i, a, ":");
        score=a[3];
        break;
      }
    }
    print read, ref, pos, score;
  }' > "$outfile"
}

# Run extraction for each file
extract_secondary_info "$SAM_BASIC" "$OUT_BASIC"
extract_secondary_info "$SAM_ALL" "$OUT_ALL"
extract_secondary_info "$SAM_IGOR" "$OUT_IGOR"

echo "✅ Secondary alignment summaries generated:"
echo "  - $OUT_BASIC"
echo "  - $OUT_ALL"
echo "  - $OUT_IGOR"

