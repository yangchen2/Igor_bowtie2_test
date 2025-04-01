# Benchmarking Bowtie2 Alignment Methods
## Objective
The primary goal is to assess the performance and biological interpretability of the genome-independent alignment method developed by [Igor via the seed_filter_large_250102 branch](https://github.com/sfiligoi/bowtie2/tree/seed_filter_large_250102). Igor's approach addresses limitations in traditional read mapping workflows arising from how Bowtie2 handles multi-mapped reads.
Traditionally, Bowtie2 captures only up to 16 multi-mapped locations per read, leading to several issues:
We do not capture all multi-mapped locations for every genome.
As reference databases grow, the likelihood of capturing all true multi-maps decreases.
Bowtie2 assumes that all genomes must be loaded into memory simultaneously, placing an artificial upper bound on the number of genomes that can be indexed.
To overcome these limitations, Igor’s genome-independent paradigm modifies Bowtie2 to: 
Capture all multi-mapped reads per genome rather than globally
It allows index partitioning across genomes, enabling infinite scaling of reference databases. 
Make it possible to retain genome-specific mappings indefinitely (e.g., for Qiita), avoiding the need to remap every sample when the reference database is updated.

## Benchmarking Strategy

This will be evaluated using real metagenomic data, focusing on the alignment behavior and resulting downstream analyses. To establish a controlled baseline, we selected a monoculture dataset containing only _Staphylococcus aureus_ (_S. aureus_) published in the [Zebra paper](https://www.sciencedirect.com/org/science/article/pii/S2379507722003592). This mock community dataset provides a reliable ground truth for assessing alignment accuracy, sensitivity, specificity, and downstream taxonomic profiling. The raw data was downloaded from [QIITA study 11919](https://qiita.ucsd.edu/study/description/11919), and only the _S. aureus_ monoculture samples will be used.

We will conduct the following alignments:

1. Baseline alignment using Bowtie2 (default version)
Input: _S. aureus_ monoculture dataset
Reference: NCBI _S. aureus_ reference genome
Purpose: Establish baseline mapping quality and rate using standard settings.


2. Genome-independent alignment using Igor’s Bowtie2 branch
Input: _S. aureus_ monoculture dataset
Reference: NCBI _S. aureus_ reference genome
Purpose: Compare accuracy and alignment rates to baseline; test genome-independent modifications.


3. Database-wide alignment using standard Bowtie2
Input: _S. aureus_ monoculture dataset
Reference: WoL2 database (comprehensive microbial reference)
Purpose: Evaluate behavior when aligning against a large, diverse database with the potential for off-target hits.


4. Database-wide alignment using Igor’s Bowtie2 branch
Input: _S. aureus_ monoculture dataset
Reference: WoL2 database
Purpose: Assess whether genome-independent alignment improves or degrades mapping accuracy, multi-mapping handling, and taxonomic resolution.
