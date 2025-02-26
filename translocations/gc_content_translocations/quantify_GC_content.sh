#!/bin/bash
#SBATCH --partition 128x24
#SBATCH --job-name=bedtools
#SBATCH --output=bedtools_%j.out
#SBATCH --error=bedtools_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kuckele@ucsc.edu
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB

# Load modules
module load samtools/1.21
module load bedtools/2.26.0

# Index FASTA file
samtools faidx C.lasius_genome_NCBI_v2.fasta

# Extract sequence data
bedtools getfasta -fi C.lasius_genome_NCBI_v2.fasta -bed roi.bed -fo roi_sequences.fasta

# Compute GC content with sliding window approach
bedtools makewindows -g C.lasius_genome_NCBI_v2.fasta.fai -w 100000 > windows.bed
bedtools nuc -fi C.lasius_genome_NCBI_v2.fasta -bed windows.bed > gc_content.txt
