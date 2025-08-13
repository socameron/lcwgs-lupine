#!/usr/bin/env Rscript
# scripts/plot_admixture.R
#
# Conventional STRUCTURE-style barplot from a PCAngsd .Q file

library(ggplot2)
library(dplyr)
library(tidyr)

# 1) parse command-line args
args     <- commandArgs(trailingOnly = TRUE)
qfile    <- args[1]  # e.g. results/pcangsd/.../all_popln_canonical_SNP_pcangsd.admix.18.Q
meta_csv <- args[2]  # data/lists/.../all_samples_to_popln_bam_order.csv
out_png  <- args[3]  # e.g. results/plots/.../all_popln_SNP_admix.png

# 2) load Q-matrix and metadata
Qmat <- read.table(qfile, header = FALSE, stringsAsFactors = FALSE)
meta <- read.csv(meta_csv, stringsAsFactors = FALSE)

# 3) annotate rows
Qmat$Individual     <- meta$Sample_name    # unique IDs, e.g. "PBBT-22"
Qmat$PopulationCode <- meta$Site_code2     # e.g. "E.MN.2"

# 4) pivot to long form
admix_long <- Qmat %>%
  mutate(Row = row_number()) %>%
  pivot_longer(
    cols      = starts_with("V"),
    names_to  = "Cluster",
    values_to = "Proportion"
  ) %>%
  select(Individual, PopulationCode, Cluster, Proportion)

# 5) order individuals by Population then by their strongest ancestry
order_df <- admix_long %>%
  group_by(Individual, PopulationCode) %>%
  summarize(MaxProp = max(Proportion), .groups="drop") %>%
  arrange(PopulationCode, desc(MaxProp))

admix_long$Individual <- factor(admix_long$Individual,
                                levels = order_df$Individual)

# 6) choose a palette for the K clusters
K <- length(unique(admix_long$Cluster))
cluster_cols <- colorRampPalette(RColorBrewer::brewer.pal(8,"Set2"))(K)

# 7) theme tweaks
custom_theme <- theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x       = element_blank(),  # hide sample labels
    axis.ticks.x      = element_blank(),
    axis.title.x      = element_blank(),
    axis.title.y      = element_text(size=12),
    plot.title        = element_text(hjust=0.5),
    legend.position   = "top"
  )

# 8) plot
p <- ggplot(admix_long, aes(x = Individual, y = Proportion, fill = Cluster)) +
  geom_col(width = 1) +
  facet_grid(~ PopulationCode, scales="free_x", space="free_x") +
  scale_fill_manual(
    values = cluster_cols,
    name   = "Cluster"
  ) +
  labs(
    title = paste0("Admixture proportions (K = ", K, ")"),
    y     = "Proportion"
  ) +
  custom_theme

# 9) save
ggsave(out_png, p, width = 14, height = 6, dpi = 300)
