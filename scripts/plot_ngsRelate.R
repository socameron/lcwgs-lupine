# Load necessary libraries
library(ggplot2)
library(dplyr)


# Define the old population IDs (from your Snakemake wildcard)
old_populations <- c("HPW", "IDNP-MW", "LCTGP", "MFNP", "PPP", "RLPLV", "SWCP", "APB", "BSL", "BSNA", "CPB", "FMB", "GRAY", "NBWA", "NGP", "PBBT", "RSB", "UWA")
populations <- old_populations

# Define a mapping from old IDs to new cluster names
pop_mapping <- c(
  "NGP"     = "C.IL.1",
  "IDNP-MW" = "C.IN.2",
  "MFNP"    = "C.MI.2",
  "RSB"     = "C.MN.4",
  "LCTGP"   = "C.OH.1",
  "SWCP"    = "C.ON.3",
  "UWA"     = "C.WI.1",
  "FMB"     = "C.WI.2",
  "GRAY"    = "E.MI.7",
  "PBBT"    = "E.MN.2",
  "BSL"     = "E.MN.3",
  "CPB"     = "E.NH.1",
  "APB"     = "E.NY.1",
  "HPW"     = "E.ON.5",
  "RLPLV"   = "E.ON.8",
  "PPP"     = "E.ON.9",
  "BSNA"    = "E.WI.4",
  "NBWA"    = "E.WI.9"
)

# Define populations and color schemes
population_colors <- c(
  # Core populations:
  "C.OH.1" = "#DA70D6",    # orchid
  "C.MI.2" = "#9932CC",    # darkorchid
  "C.ON.3" = "#551A8B",    # purple4
  "C.IN.2" = "#B19CD9",    # softer purple
  "C.IL.1" = "#9370DB",    # mediumpurple
  "C.MN.4" = "#BA55D3",    # mediumorchid
  "C.WI.1" = "#663399",    # rebeccapurple
  "C.WI.2" = "#9400D3",    # darkviolet
  
  #Edge populations
  "E.ON.5" = "#00FFFF",    # cyan
  "E.ON.8" = "#7FFFD4",    # aquamarine
  "E.ON.9" = "#00868B",    # turquoise4
  "E.MI.7" = "#20B2AA",    # lightseagreen
  "E.MN.2" = "#48D1CC",    # mediumturquoise
  "E.MN.3" = "#5F9EA0",    # cadetblue
  "E.NH.1" = "#008B8B",    # darkcyan
  "E.NY.1" = "#00CED1",    # dark turqoise 
  "E.WI.4" = "#40E0D0",    # turquoise
  "E.WI.9" = "#009999")     # balanced teal

# Create popln labels from color code
new_populations <- unnames(pop_mapping[old_populations])

# Initialize empty data frame
combined_theta_data <- data.frame()
# Initialize an empty data frame to store the combined rab values from all populations
combined_data <- data.frame()

# Loop through populations and read the corresponding ngsRelate output files
for (pop_old in populations) {
  # Construct file path
  file_path <- paste0("results/ngsRelate/hap2/", pop_old, "_ngsrelate.out")
  
  # Read the ngsRelate output file
  ngsrelate_data <- read.table(file_path, header = TRUE)
  
  # Convert old popln names to new
  new_pop <- pop_mapping[pop_old]

  # Replace population names
  pop_group <- ifelse(startsWith(new_pop, "C."), "Core", "Edge")
  
  # Extract only the 'rab' column and add a population identifier
  pop_data <- ngsrelate_data %>%
    select(rab) %>%
    mutate(population = new_pop,
           population_group = pop_group)
  
  # Combine the population data into the main data frame
  combined_data <- rbind(combined_data, pop_data)
}

# Reorder populations by core and edge
combined_data$population <- factor(combined_data$population, levels = c(
  "C.IN.2", "C.OH.1", "C.MI.2", "C.ON.3", "C.IL.1", "C.MN.4", "C.WI.1", "C.WI.2",
  "E.ON.5", "E.ON.8", "E.ON.9", "E.MI.7", "E.MN.2", "E.MN.3", "E.NH.1", "E.NY.1", "E.WI.4", "E.WI.9"
))

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
