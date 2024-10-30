# Loading necessary libraries
library(ggplot2)

# Process LR files
args <- commandArgs(trailingOnly = TRUE)
print(args)
input_file <- args[1]
pop_names <- args[2]
plot_file <- args[3]


input_file <- "results/pcangsd/hap2/canonical/pcangsd_input/all_popln_canonical_SNP_pcangsd.cov"
pop_names <- "data/lists/Batch1_popln_to_samples.info"
plot_file <- "results/plots/hap2/PCAngsd/all_popln_SNP_PCA.png"


# Read population names and covariance matrix
pop <- read.table(pop_names, header = FALSE, sep = "\t")
C <- as.matrix(read.table(input_file, header = FALSE))
e <- eigen(C)


# Prepare data for ggplot
pca_data <- data.frame(PC1 = e$vectors[, 1], PC2 = e$vectors[, 2], Population = pop[, 1])

# Define colors for each population-based cluster (same as in admixture plot)
cluster_colors <- c("C.OH.1" = "orchid", "C.MI.2" = "darkorchid", 
                    "C.ON.3" = "purple4", "C.IN.2" = "lightsalmon1", 
                    "E.ON.5" = "cyan", "E.ON.8" = "aquamarine", 
                    "E.ON.9" = "turquoise4")

# Set the Population column as a factor with specific levels to maintain consistent coloring
pca_data$Population <- factor(pca_data$Population, levels = names(cluster_colors))

# Define a custom theme
custom_theme <- theme(panel.background = element_rect(fill = "white", color = NA),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black", linewidth = 1), 
                      axis.ticks = element_line(linewidth = 1), 
                      axis.ticks.length = unit(0.25, "cm"), 
                      plot.title = element_text(hjust = 0.5), 
                      axis.text = element_text(size=12),
                      axis.title = element_text(size=12),
                      legend.position = "top", 
                      legend.margin = margin(0, 0, 0, 0), 
                      legend.text = element_text(size=12), 
                      legend.title = element_text(size=10))

# Create the plot using ggplot
gg_pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(size = 3) +
  labs(title = "Individual Allele Frequency", x = "PC1", y = "PC2") +
  scale_color_manual(values = cluster_colors, name = "Population") +
  custom_theme

# Save the plot as a PNG file
ggsave(plot_file, plot = gg_pca_plot, width = 10, height = 8)
