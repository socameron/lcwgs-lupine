#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 3) {
  stop("Usage: Rscript flag_admixed_individuals.R <qopt_file> <bamlist_file> <metadata_csv>")
}

qopt_fp <- args[1]
bamlist_fp <- args[2]
meta_fp <- args[3]

# Load input
qmat <- read.table(qopt_fp, header = FALSE)
bamlist <- readLines(bamlist_fp)
meta <- read.csv(meta_fp)

# Extract sample names from BAM paths
sample_names <- sub("_hap2_realign.bam$", "", basename(bamlist))

# Combine
if(length(sample_names) != nrow(qmat)) {
  stop("Mismatch: number of samples in BAM list and .qopt file")
}
qmat$Sample_name <- sample_names

# Merge with metadata
df <- merge(qmat, meta, by = "Sample_name", all.x = TRUE)

# Calculate max Q per individual
df$maxQ <- apply(df[, grep("^V", names(df))], 1, max)

# Flag admixed
df$admixed <- df$maxQ < 0.9

# Split by population
library(dplyr)
summary_df <- df %>%
  group_by(Site_code) %>%
  arrange(desc(maxQ)) %>%
  select(Sample_name, Site_code, maxQ, admixed)

# Write report
write.table(summary_df, file = "results/NeEstimator/admixture_q90_flagged.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done: Admixture report written to results/NeEstimator/admixture_q90_flagged.txt\n")
