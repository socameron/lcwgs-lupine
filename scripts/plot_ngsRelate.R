# Load necessary libraries
library(ggplot2)
library(dplyr)

# Replace old population names with new ones
old_levels <- c("ARCA", "RLPLV", "PPP", "HPW", "LW-KBS", "IDNP-MW", "LCTGP", "MFNP", "NCBW-FOT", "SWCP")
new_levels <- c("E.MI.3", "E.ON.8", "E.ON.9", "E.ON.5", "E.ON.12", "C.IN.2", "C.OH.1", "C.MI.2", "C.ON.1", "C.ON.3")

# Define populations and color scheme
core_populations <- c("LCTGP", "MFNP", "IDNP-MW", "SWCP")
edge_populations <- c("PPP", "HPW", "RLPLV")
populations <- c(core_populations, edge_populations)

# Create a color mapping for core and edge populations
population_colors <- c("darkorchid3", "darkorchid3", "darkorchid3", "darkorchid3", # Core colors
                       "cyan3", "cyan3", "cyan3") # Edge colors

# Initialize an empty data frame to store the combined rab values from all populations
combined_data <- data.frame()

# Loop through populations and read the corresponding ngsRelate output files
for (pop in populations) {
  # Construct file path
  file_path <- paste0("results/ngsRelate/hap2/", pop, "_ngsrelate.out")
  
  # Read the ngsRelate output file
  ngsrelate_data <- read.table(file_path, header = TRUE)
  
  # Replace population names
  pop <- sapply(pop, function(x) {
    idx <- match(x, old_levels)
    if (!is.na(idx)) new_levels[idx] else x
  })
  
  # Extract only the 'rab' column and add a population identifier
  pop_data <- ngsrelate_data %>%
    select(rab) %>% # Remove rows where rab is 0
    mutate(population = pop,
           population_group = ifelse(pop %in% core_populations, "Core", "Edge"))
  
  # Combine the population data into the main data frame
  combined_data <- rbind(combined_data, pop_data)
}

# Reorder populations by core and edge
combined_data$population <- factor(combined_data$population, levels = c("C.IN.2", "C.OH.1", "C.MI.2", "C.ON.1", "C.ON.3", "E.MI.3", "E.ON.5", "E.ON.8", "E.ON.9", "E.ON.12"))

# Custom theme settings
custom_theme <- theme(panel.border = element_blank(), 
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black", linewidth = 1), 
                      axis.ticks = element_line(linewidth = 1), 
                      axis.ticks.length = unit(0.25, "cm"), 
                      plot.title = element_text(hjust = 0.5), 
                      axis.text = element_text(size=12, angle = 45, hjust = 1), 
                      axis.title = element_text(size=12),
                      legend.position = "none", 
                      legend.margin = margin(0, 0, 0, 0), 
                      legend.text = element_text(size=12), 
                      legend.title = element_text(size=10))

# Create dataset for non-zeroes
non_zero_data <- combined_data %>% filter(rab > 0)

# Boxplot with zeroes
boxplot_1 <- ggplot(combined_data, aes(x = population, y = rab, fill = population)) +
  geom_boxplot(aes(alpha=0.3)) +
  scale_fill_manual(values = population_colors) +
  labs(title = "Boxplot of Relatedness (rab) Across Populations",
       x = "Population",
       y = "Pairwise Relatedness") +
  custom_theme
ggsave("results/plots/hap2/ngsRelate/rab_boxplot.png", plot = boxplot_1, width = 8, height = 6)

# Boxplot without zeroes
boxplot_2 <- ggplot(non_zero_data, aes(x = population, y = rab, fill = population)) +
  geom_boxplot(aes(alpha=0.3)) +
  scale_fill_manual(values = population_colors) +
  labs(title = "Boxplot of Relatedness (rab) Across Populations",
       x = "Population",
       y = "Pairwise Relatedness") +
  custom_theme
ggsave("results/plots/hap2/ngsRelate/rab_boxplot_no_zeros.png", plot = boxplot_2, width = 8, height = 6)

# Boxplot without zeroes with log transform
boxplot_3 <- ggplot(non_zero_data, aes(x = population, y = log10(rab), fill = population)) +
  geom_boxplot(aes(alpha=0.3)) +
  scale_fill_manual(values = population_colors) +
  labs(title = "Boxplot of Relatedness (rab) Across Populations",
       x = "Population",
       y = "log10 Pairwise Relatedness") +
  custom_theme
ggsave("results/plots/hap2/ngsRelate/rab_log_boxplot_no_zeros.png", plot = boxplot_3, width = 8, height = 6)

# Violin Plot
violin_plot <- ggplot(non_zero_data, aes(x = population, y = rab, fill = population)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = population_colors) +
  labs(title = "Violin Plot of Relatedness (rab) Across Populations",
       x = "Population",
       y = "Pairwise Relatedness") +
  custom_theme
ggsave("results/plots/hap2/ngsRelate/rab_violin_plot_no_zeros.png", plot = violin_plot, width = 8, height = 6)

# Split Violin Plot with Boxplot Overlay
split_violin_plot <- ggplot(non_zero_data, aes(x = population, y = rab, fill = population_group)) +
  geom_violin(trim = FALSE, aes(fill = population_group), position = position_dodge(0.9)) + # Split violin plot
  geom_boxplot(width = 0.2, outlier.shape = NA, position = position_dodge(0.9), color = "black") + # Overlay the boxplot
  scale_fill_manual(values = c("Core" = "darkorchid3", "Edge" = "cyan3")) + # Adjust colors for Core and Edge
  labs(title = "Split Violin and Boxplot of Relatedness (rab) Across Populations (Excluding Zeros)",
       x = "Population",
       y = "Pairwise Relatedness") +
  custom_theme
ggsave("results/plots/hap2/ngsRelate/rab_split_violin_boxplot_no_zeros.png", plot = split_violin_plot, width = 8, height = 6)
