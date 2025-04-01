import os
import pysam
import pandas as pd
import matplotlib.pyplot as plt

# Define directories
basic_dir = "output_sam_files_basic"
igor_dir = "output_sam_files_igor_fixed"
output_dir = "slurm_out"
os.makedirs(output_dir, exist_ok=True)

# Function to compute alignment statistics
def parse_sam(file):
    total_reads = 0
    mapped_reads = 0
    properly_paired = 0
    unique_mapped = 0
    mapq_scores = []

    try:
        with pysam.AlignmentFile(file, "r", check_sq=False) as samfile:  # <-- Added check_sq=False
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

    except ValueError as e:
        print(f"Error reading {file}: {e}")
        return {
            "Total Reads": 0,
            "Mapped Reads": 0,
            "Properly Paired Reads": 0,
            "Uniquely Mapped Reads": 0,
            "Average MAPQ": 0
        }

#list of SAM files in both directories
basic_sam_files = {f.split(".sam")[0]: os.path.join(basic_dir, f) for f in os.listdir(basic_dir) if f.endswith(".sam")}
igor_sam_files = {f.split("_l100_d.sam")[0]: os.path.join(igor_dir, f) for f in os.listdir(igor_dir) if f.endswith(".sam")}

# Find matching samples
common_samples = set(basic_sam_files.keys()) & set(igor_sam_files.keys())

# Process each sample and compare statistics
stats_list = []
for sample in common_samples:
    basic_stats = parse_sam(basic_sam_files[sample])
    igor_stats = parse_sam(igor_sam_files[sample])

    stats_list.append({"Sample": sample, "Method": "Regular Bowtie2", **basic_stats})
    stats_list.append({"Sample": sample, "Method": "Igor's Bowtie2", **igor_stats})

# Convert to DataFrame
df = pd.DataFrame(stats_list)

# Save results to file
comparison_file = os.path.join(output_dir, "alignment_comparison.csv")
df.to_csv(comparison_file, index=False)

# Print summary
print(df)

# Plot comparison
metrics = ["Mapped Reads", "Properly Paired Reads", "Uniquely Mapped Reads", "Average MAPQ"]
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

for ax, metric in zip(axes.flatten(), metrics):
    df.pivot(index="Sample", columns="Method", values=metric).plot(kind="bar", ax=ax, alpha=0.75)
    ax.set_title(metric)
    ax.set_ylabel("Count" if "Reads" in metric else "Score")
    ax.grid(axis="y", linestyle="--", alpha=0.7)

plt.suptitle("Comparison of Bowtie2 Alignment Methods")
plt.tight_layout()

# Save the figure
plot_file = os.path.join(output_dir, "alignment_comparison.png")
plt.savefig(plot_file, dpi=300)

print(f"Results saved to {comparison_file}")
print(f"Plot saved to {plot_file}")

