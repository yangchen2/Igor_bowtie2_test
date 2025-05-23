import os
import pysam
import pandas as pd
import matplotlib.pyplot as plt

# Define input directories
# basic_dir = "output_sam_files_basic"
# basic_all_reads_dir = "output_sam_files_basic_all-reads"
# igor_dir = "output_sam_files_igor"
# output_dir = "comparison_results"

basic_dir = "output_sam_files_basic_wolr2"
basic_all_reads_dir = "output_sam_files_basic_all-reads_wolr2"
igor_dir = "output_sam_files_igor_wolr2"
output_dir = "comparison_results/monoculture_wolr2"
os.makedirs(output_dir, exist_ok=True)

# Function to compute alignment statistics
def parse_sam(file):
    total_reads = 0
    mapped_reads = 0
    properly_paired = 0
    unique_mapped = 0
    mapq_scores = []

    try:
        with pysam.AlignmentFile(file, "r", check_sq=False) as samfile:
            for read in samfile:
                total_reads += 1
                if not read.is_unmapped:
                    mapped_reads += 1
                    mapq_scores.append(read.mapping_quality)
                    if read.mapping_quality > 0:
                        unique_mapped += 1
                    if read.is_proper_pair:
                        properly_paired += 1

        avg_mapq = sum(mapq_scores) / len(mapq_scores) if mapq_scores else 0

        return {
            "Total Reads": total_reads,
            "Mapped Reads": mapped_reads,
            "Properly Paired Reads": properly_paired,
            "Uniquely Mapped Reads": unique_mapped,
            "Average MAPQ": avg_mapq
        }

    except Exception as e:
        print(f"Error reading {file}: {e}")
        return {
            "Total Reads": 0,
            "Mapped Reads": 0,
            "Properly Paired Reads": 0,
            "Uniquely Mapped Reads": 0,
            "Average MAPQ": 0
        }

# Collect files for all three methods
basic_sam_files = {f.split(".sam")[0]: os.path.join(basic_dir, f) for f in os.listdir(basic_dir) if f.endswith(".sam")}
basic_all_reads_sam_files = {f.split(".sam")[0]: os.path.join(basic_all_reads_dir, f) for f in os.listdir(basic_all_reads_dir) if f.endswith(".sam")}
igor_sam_files = {f.split("_l100_d.sam")[0]: os.path.join(igor_dir, f) for f in os.listdir(igor_dir) if f.endswith(".sam")}

# Find matching sample names
common_samples = set(basic_sam_files) & set(basic_all_reads_sam_files) & set(igor_sam_files)

# Parse and collect statistics
stats_list = []
for sample in sorted(common_samples):
    stats_list.append({
        "Sample": sample,
        "Method": "Standard Bowtie2 (-k 16)",
        **parse_sam(basic_sam_files[sample])
    })
    stats_list.append({
        "Sample": sample,
        "Method": "Standard Bowtie2 (-a)",
        **parse_sam(basic_all_reads_sam_files[sample])
    })
    stats_list.append({
        "Sample": sample,
        "Method": "Igor's Bowtie2 (-a -d -l 100)",
        **parse_sam(igor_sam_files[sample])
    })

# Convert to DataFrame
df = pd.DataFrame(stats_list)

# Save to CSV
comparison_file = os.path.join(output_dir, "alignment_comparison_all_method_wolr2.csv")
df.to_csv(comparison_file, index=False)
print(f"Saved comparison table to {comparison_file}")

# Plotting
metrics = ["Mapped Reads", "Properly Paired Reads", "Uniquely Mapped Reads", "Average MAPQ"]
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

for ax, metric in zip(axes.flatten(), metrics):
    df.pivot(index="Sample", columns="Method", values=metric).plot(kind="bar", ax=ax, alpha=0.8)
    ax.set_title(metric)
    ax.set_ylabel("Count" if "Reads" in metric else "Score")
    ax.set_xlabel("Sample")
    ax.grid(axis="y", linestyle="--", alpha=0.6)

plt.suptitle("Comparison of Bowtie2 Alignment Methods (All Variants) Against Wolr2", fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plot_file = os.path.join(output_dir, "alignment_comparison_all_methods_wolr2.png")
plt.savefig(plot_file, dpi=300)
print(f"Saved plot to {plot_file}")
