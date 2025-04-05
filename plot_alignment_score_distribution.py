import os
import pysam
import matplotlib.pyplot as plt

# BAM paths
# bam_standard = "output_sam_files_basic/bams/NC.Staphy2.1A_sorted.bam"
bam_standard = "output_sam_files_basic_all-reads/bams/NC.Staphy2.1A_sorted.bam"

bam_igor = "output_sam_files_igor_fixed/bams/NC.Staphy2.1A_fixed_sorted.bam"

def extract_alignment_scores(bam_path):
    scores = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if not read.is_unmapped:
                for tag in read.get_tags():
                    if tag[0] == "AS":
                        scores.append(tag[1])
                        break
    return scores

scores_standard = extract_alignment_scores(bam_standard)
scores_igor = extract_alignment_scores(bam_igor)

# Plot
plt.figure(figsize=(10, 5))
plt.hist(scores_standard, bins=30, alpha=0.6, label="Standard Bowtie2 with -a flag", edgecolor='black')
plt.hist(scores_igor, bins=30, alpha=0.6, label="Igor Bowtie2", edgecolor='black')
plt.title("Alignment Score Distribution: NC.Staphy2.1A")
plt.xlabel("Alignment Score (AS:i)")
plt.ylabel("Read Count")
plt.legend()
plt.tight_layout()
plt.savefig("alignment_score_comparison_NC.Staphy2.1A_all-standard.png")