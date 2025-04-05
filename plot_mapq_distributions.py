import os
import pysam
import matplotlib.pyplot as plt
from collections import defaultdict

# Set input and output directories
bam_dir = "output_sam_files_igor/bams"  # your directory of BAM files
output_dir = "output_sam_files_igor/bams/mapq_plots"
os.makedirs(output_dir, exist_ok=True)

# Loop through all BAM files
for filename in os.listdir(bam_dir):
    if filename.endswith(".bam"):
        bam_path = os.path.join(bam_dir, filename)
        sample_name = os.path.splitext(filename)[0]
        mapq_counts = defaultdict(int)

        # Read BAM and count MAPQ scores
        with pysam.AlignmentFile(bam_path, "rb") as bamfile:
            for read in bamfile:
                if not read.is_unmapped:
                    mapq_counts[read.mapping_quality] += 1

        # Sort MAPQ values
        sorted_mapq = sorted(mapq_counts.items())

        # Prepare data for plotting
        x = [k for k, _ in sorted_mapq]
        y = [v for _, v in sorted_mapq]

        # Plot
        plt.figure(figsize=(8, 5))
        plt.bar(x, y, width=2, edgecolor='black')
        plt.title(f'MAPQ Distribution: {sample_name}')
        plt.xlabel('MAPQ Score')
        plt.ylabel('Read Count')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{sample_name}_MAPQ.png"))
        plt.close()

        print(f"Saved plot: {sample_name}_MAPQ.png")
