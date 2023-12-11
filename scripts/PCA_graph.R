# R script to create a graph from "all_popln_canonical_SNP_PCAngsd.cov" data

# Process LR files
args <- commandArgs(trailingOnly = TRUE)
print(args)
input_file <- args[1]
pop_names <- args[2]
plot_file <- args[3]

# Read population names and covariance matrix
pop <- read.table(pop_names, header = FALSE, sep = "\t")
C <- as.matrix(read.table(input_file, header = FALSE))
e <- eigen(C)

# Loading ggplot for plot production
library(ggplot2)

# Prepare data for ggplot
pca_data <- data.frame(PC1 = e$vectors[, 1], PC2 = e$vectors[, 2], Population = pop[, 1])

# Create the plot using ggplot
gg_pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Population)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Individual Allele Frequency", x = "PC1", y = "PC2") +
  scale_color_discrete(name = "Population")

# Save the plot as a PNG file
ggsave(plot_file, plot = gg_pca_plot, width = 10, height = 8)