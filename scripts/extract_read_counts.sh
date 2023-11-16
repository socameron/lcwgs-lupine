#!/bin/bash

# Redirect stdout and stderr to log file
exec &> >(tee -a "results/logs/extract_read_counts.log")

# Output file
output_file="results/qc/counts/read_counts.txt"

# Header
echo -e "File Name\tSet Type\tNumber of Reads" > "$output_file"

# Iterate over prefilter BAM files
for bam_file in results/bam_mkdup/hap1/*_hap1_mkdup.bam; do
    base_name=$(basename "$bam_file")
    echo "Processing $base_name"
    read_count=$(module load samtools/1.17 && samtools view -c "$bam_file" -@ 2)
    echo "Read count for $base_name: $read_count"
    echo -e "$base_name\tbam_mkdup\t$read_count" >> "$output_file"
done

# Iterate over postfilter BAM files
for bam_file in results/bam_ngs_prep/hap1/*_ngs_prep.bam; do
    base_name=$(basename "$bam_file")
    echo "Processing $base_name"
    read_count=$(module load samtools/1.17 && samtools view -c "$bam_file" -@ 2)
    echo "Read count for $base_name: $read_count"
    echo -e "$base_name\tbam_ngs_prep\t$read_count" >> "$output_file"
done

