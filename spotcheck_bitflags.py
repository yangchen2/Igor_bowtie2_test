import os
import pysam
from collections import Counter, defaultdict

# bam_path = "output_sam_files_igor_fixed/bams/NC.Staphy2.1A_fixed_sorted.bam"  # update as needed
# bam_path = "output_sam_files_basic/bams/NC.Staphy2.1A_sorted.bam"  # update as needed
bam_path = "output_sam_files_basic_all-reads/bams/NC.Staphy2.1A_sorted.bam"  # update as needed

num_reads_per_flag = 3

flag_counts = Counter()
example_reads = defaultdict(list)

with pysam.AlignmentFile(bam_path, "rb") as bam:
    for read in bam:
        flag = read.flag
        flag_counts[flag] += 1
        if len(example_reads[flag]) < num_reads_per_flag:
            example_reads[flag].append(read)

print(f"\nTop 10 most common bitflags in {os.path.basename(bam_path)}:")
for flag, count in flag_counts.most_common(10):
    print(f"FLAG: {flag:<5} | Count: {count}")

print("\n--- Example reads per flag ---")
for flag, reads in list(example_reads.items())[:5]:  # limit to top 5 for display
    print(f"\nFLAG {flag}:")
    for read in reads:
        print(f"Read name: {read.query_name}, is_unmapped: {read.is_unmapped}, is_secondary: {read.is_secondary}, is_supplementary: {read.is_supplementary}")
