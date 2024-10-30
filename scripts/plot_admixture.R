# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
info_file <- args[2]
output_file <- args[3]

# Read the .Q file into a data frame
admix_data <- read.table(input_file, header = FALSE)

# Read the population info file
pop_info <- read.table(info_file, header = FALSE, stringsAsFactors = FALSE)
colnames(pop_info) <- c("PopulationCode", "Individual")

# Add individual IDs based on the order (assuming the same order as in pop_info)
admix_data <- admix_data %>%
  mutate(Individual = pop_info$Individual,
         PopulationCode = pop_info$PopulationCode)

# Calculate the average proportion per cluster for each population
cluster_means <- admix_data %>%
  pivot_longer(cols = starts_with("V"), names_to = "Cluster", values_to = "Proportion") %>%
  group_by(PopulationCode, Cluster) %>%
  summarize(AverageProportion = mean(Proportion), .groups = "drop")

# Identify the dominant cluster for each population
dominant_clusters <- cluster_means %>%
  group_by(PopulationCode) %>%
  slice_max(order_by = AverageProportion, n = 1) %>%
  select(PopulationCode, DominantCluster = Cluster)

# Create a mapping between Cluster (V1, V2, etc.) and PopulationCode based on the dominant cluster
cluster_map <- setNames(dominant_clusters$PopulationCode, dominant_clusters$DominantCluster)

# Join the dominant cluster information back to the admixture data
admix_data <- admix_data %>%
  left_join(dominant_clusters, by = "PopulationCode")

# Convert the data frame into long format for ggplot2
admix_long <- admix_data %>%
  pivot_longer(cols = starts_with("V"), names_to = "Cluster", values_to = "Proportion")

# Map the generic cluster names (V1, V2, ...) to population-based labels using the cluster_map
admix_long$Cluster <- factor(admix_long$Cluster, levels = names(cluster_map), labels = cluster_map)

# Define colors for each population-based cluster (adjust according to your populations)
cluster_colors <- c("C.OH.1" = "orchid", "C.MI.2" = "darkorchid", 
                    "C.ON.3" = "purple4", "C.IN.2" = "lightsalmon1", 
                    "E.ON.5" = "cyan", "E.ON.8" = "aquamarine", 
                    "E.ON.9" = "turquoise4")

# Set the order of the Individuals by PopulationCode and within each population by the dominant proportion
admix_long <- admix_long %>%
  group_by(Individual, PopulationCode) %>%
  mutate(DominantProportion = max(Proportion)) %>%
  ungroup()

# Order the individuals: first by PopulationCode and then by DominantProportion within each population
admix_long <- admix_long %>%
  arrange(PopulationCode, desc(DominantProportion))

# Define a custom theme
custom_theme <- theme(panel.background = element_rect(fill = "white", color = NA),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black", linewidth = 1), 
                      axis.ticks = element_line(linewidth = 1), 
                      axis.ticks.length = unit(0.25, "cm"), 
                      plot.title = element_text(hjust = 0.5), 
                      axis.text.y = element_text(size=12),
                      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
                      axis.title = element_text(size=12),
                      legend.position = "top", 
                      legend.margin = margin(0, 0, 0, 0), 
                      legend.text = element_text(size=12), 
                      legend.title = element_text(size=10))

# Number of clusters (K) is the number of columns in admix_data excluding Individual and PopulationCode
K <- ncol(admix_data) - 3

# Create the admixture bar plot with colors based on the population-aligned clusters
p <- ggplot(admix_long, aes(x = Individual, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack", width = 1) +
  scale_fill_manual(values = cluster_colors) +
  labs(x = "Individuals", y = "Admixture Proportion", title = paste("Admixture Proportions (K =", K, ")")) +
  custom_theme

# Save the plot to the specified output file
ggsave(output_file, plot = p, width = 14, height = 8)
