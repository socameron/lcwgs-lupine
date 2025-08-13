#!/usr/bin/env Rscript
# Usage: Rscript plot_ld_distribution.R popA_ld_sample.txt popB_ld_sample.txt ... out.pdf

args <- commandArgs(trailingOnly=TRUE)
out_plot <- tail(args, n=1)
in_files <- head(args, -1)

library(ggplot2)
library(dplyr)
library(readr)

# Read and tag each sample file
ld_list <- lapply(in_files, function(f){
  pop <- sub("_ld_sample.txt$", "", basename(f))
  # read_tsv will read a single-column file into a tibble
  tibble(
    population = pop,
    r2 = read_tsv(f, col_names="r2", show_col_types=FALSE)$r2
  )
})

# Combine into one data frame
ld_df <- bind_rows(ld_list)

# Plot: histograms or density
p <- ggplot(ld_df, aes(x = r2, fill = population)) +
  geom_histogram(bins = 50, alpha=0.6, position="identity") +
  facet_wrap(~population, scales="free_y", ncol=4) +
  labs(x = expression(r^2), y = "Count",
       title = "LD (r²) Distribution by Population (sampled)") +
  theme_minimal() +
  theme(legend.position="none")

ggsave(out_plot, p, width=12, height=8)
