# Load necessary libraries
library(ggplot2)
library(dplyr)


# Define the old population IDs (from your Snakemake wildcard)
old_populations <- c("HPW", "IDNP-MW", "LCTGP", "MFNP", "PPP", "RLPLV", "SWCP", 
                    "APB", "BSL", "BSNA", "CPB", "FMB", "GRAY", "NBWA", "NGP", "PBBT", "RSB", "UWA")
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
new_populations <- unname(pop_mapping[old_populations])

# Initialize empty data frame
combined_theta_data <- data.frame()

# Function to read theta data from .pestPG files
read_theta_data <- function(file_path) {
  # Read the tab-separated file, skipping the first comment line
  theta_data <- read.table(file_path, header = TRUE, sep = "\t", comment.char = "#", skip = 1)

  # Rename columns
  colnames(theta_data) <- c("Region", "Chr", "WinCenter", "tW", "tP", "tF", "tH", "tL", 
                            "TajimaD", "FuLiF", "FuLiD", "FayH", "Zeng", "nSites")
  
  # Select relevant columns (window position, Tajima's D, and nucleotide diversity)
  theta_data <- theta_data %>% 
    select(window_position = WinCenter, tajimasD = TajimaD, nucleotide_diversity = tP, sites = nSites)
  
  # Calculate the 15th percentile cutoff for nSites within this population
  nSites_threshold <- quantile(theta_data$sites, 0.15, na.rm = TRUE)
  
  # Filter out the lowest 15% of nSites values
  theta_data <- theta_data %>% filter(sites > nSites_threshold)

  return(theta_data)
}



# Loop through populations and read their respective theta files
for (pop_old in populations) {
  file_path <- paste0("results/theta/hap2/", pop_old, "_out.thetasWindow.gz.pestPG")
  
  # Read the theta data
  theta_data <- read_theta_data(file_path)
  
  # Convert old popln names to new name
  new_pop <- pop_mapping[pop_old]

  # Replace population names
  pop_group <- ifelse(startsWith(new_pop, "C."), "Core", "Edge")
  
  # Add population identifier and group (Core/Edge)
  theta_data <- theta_data %>%
    mutate(population = new_pop,
          population_group = pop_group)
  
  # Combine data for all populations
  combined_theta_data <- rbind(combined_theta_data, theta_data)
}

# Reorder populations by core and edge
combined_theta_data$population <- factor(combined_theta_data$population, levels = c(
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
                      axis.text = element_text(size=14, angle = 45, hjust = 1), 
                      axis.title = element_text(size=14),
                      legend.position = "none", 
                      legend.margin = margin(0, 0, 0, 0), 
                      legend.text = element_text(size=14), 
                      legend.title = element_text(size=10))

mean_values <- combined_theta_data %>%
  group_by(population) %>%
  summarise(mean_diversity = mean(nucleotide_diversity/sites, na.rm = TRUE))

median_values <- combined_theta_data %>%
  group_by(population) %>%
  summarise(median_diversity = median(nucleotide_diversity/sites, na.rm = TRUE))


## BOX PLOTS

# Plotting boxplot for both Tajima's D and Nucleotide Diversity
p_tajimasD <- combined_theta_data %>%
  filter(!sites == 0) %>%
  ggplot(aes(x = population, y = tajimasD/sites, fill = population)) +
  geom_boxplot() +
  scale_fill_manual(values = population_colors) +
  scale_x_discrete(labels = function(x) sub("^[CE]\\.", "", x)) + # remove C. and E. precursors of population names
  labs(title = "Tajima's D across Populations", x = "Population", y = "Tajima's D") +
  custom_theme

p_nucleotide_div <- combined_theta_data %>%
  filter(!sites == 0) %>% 
  ggplot(aes(x = population, y = nucleotide_diversity/sites, fill = population)) +
  geom_boxplot() +
  #geom_hline(data = mean_values, aes(yintercept = mean_diversity), color = "blue", linetype = "dashed", size = 0.5) +
  #geom_hline(data = median_values, aes(yintercept = median_diversity), color = "red", linetype = "solid", size = 0.5) +
  scale_fill_manual(values = population_colors) +
  scale_x_discrete(labels = function(x) sub("^[CE]\\.", "", x)) + # remove C. and E. precursors of population names
  labs(title = "Nucleotide Diversity (π) across Populations", x = "Population", y = "Nucleotide Diversity (π)") +
  custom_theme

# Save the plots
ggsave("results/plots/hap2/theta/tajimasD_boxplot.tiff", plot = p_tajimasD, width = 10, height = 6)
ggsave("results/plots/hap2/theta/nucleotide_diversity_boxplot.tiff", plot = p_nucleotide_div, width = 10, height = 6)



## SUMMARY POINTS

# Calculate summary statistics: mean, standard deviation, and standard error for nucleotide diversity
nucleotide_summary <- combined_theta_data %>%
  filter(sites != 0) %>%  # Ensure you exclude windows with zero sites
  group_by(population) %>%
  summarise(
    mean_diversity = mean(nucleotide_diversity/sites, na.rm = TRUE),
    sd_diversity = sd(nucleotide_diversity/sites, na.rm = TRUE),
    n = n(),
    se_diversity = sd_diversity / sqrt(n)
  )

tajima_summary <- combined_theta_data %>%
  filter(sites != 0) %>%  # Ensure you exclude windows with zero sites
  group_by(population) %>%
  summarise(
    mean_tajimasD = mean(tajimasD/sites, na.rm = TRUE),
    sd_tajimasD = sd(tajimasD/sites, na.rm = TRUE),
    n = n(),
    se_tajimasD = sd_tajimasD / sqrt(n))

# Preview the summary table (optional)
print(nucleotide_summary)
print(tajima_summary)

# Create a bar plot of the mean nucleotide diversity with error bars representing the standard error
p_mean_nucleotide <- ggplot(nucleotide_summary, aes(x = population, y = mean_diversity, color = population)) +
  geom_point(size = 6) +
  geom_errorbar(aes(ymin = mean_diversity - se_diversity, ymax = mean_diversity + se_diversity), 
                width = 0.75, linewidth = 1.5) +
  scale_color_manual(values = population_colors) +
  labs(
    title = "Nucleotide Diversity at 50KB sliding window",
    x = "Population",
    y = expression("Nucleotide Diversity (π)")
  ) +
  scale_x_discrete(labels = function(x) sub("^[CE]\\.", "", x)) + # remove C. and E. precursors of population names
  custom_theme

# Create a bar plot of the mean nucleotide diversity with error bars representing the standard error
p_mean_tajima <- ggplot(tajima_summary, aes(x = population, y = mean_tajimasD, color = population)) +
  geom_point(size = 6) +
  geom_errorbar(aes(ymin = mean_tajimasD - se_tajimasD, ymax = mean_tajimasD + se_tajimasD), 
                width = 0.75, linewidth = 1.5) +
  scale_color_manual(values = population_colors) +
  labs(
    title = "Tajima's D at 50KB sliding window",
    x = "Population",
    y = expression("Tajima's D")
  ) +
  scale_x_discrete(labels = function(x) sub("^[CE]\\.", "", x)) + # remove C. and E. precursors of population names
  custom_theme

# Save the plot to file
ggsave("results/plots/hap2/theta/nucleotide_diversity_mean_se.tiff", plot = p_mean_nucleotide, width = 8, height = 6)
ggsave("results/plots/hap2/theta/tajimasD_mean_se.tiff", plot = p_mean_tajima, width = 8, height = 6)







## VIOLIN PLOTS

# Calculate the 90th percentile of the nucleotide diversity
threshold <- quantile(combined_theta_data$nucleotide_diversity, 0.90, na.rm = TRUE)

# Plotting violin plot for both Tajima's D and Nucleotide Diversity
p_tajimasD_violin <- combined_theta_data %>%
  filter(!sites == 0) %>%
  ggplot(aes(x = population, y = tajimasD/sites, fill = population)) +
  geom_violin() + 
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") +
  scale_fill_manual(values = population_colors) +
  scale_x_discrete(labels = function(x) sub("^[CE]\\.", "", x)) + # remove C. and E. precursors of population names
  labs(title = "Tajima's D across Populations", x = "Population", y = "Tajima's D") +
  custom_theme

p_nucleotide_div_violin <- combined_theta_data %>%
  filter(!sites == 0) %>%
  ggplot(aes(x = population, y = nucleotide_diversity/sites, fill = population)) +
  geom_violin(width=1) + 
  #geom_hline(data = mean_values, aes(yintercept = mean_diversity), color = "blue", linetype = "dashed", size = 0.5) +
  #geom_hline(data = median_values, aes(yintercept = median_diversity), color = "red", linetype = "solid", size = 0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") +
  scale_fill_manual(values = population_colors) +
  scale_x_discrete(labels = function(x) sub("^[CE]\\.", "", x)) + # remove C. and E. precursors of population names
  #scale_y_continuous(breaks = seq(floor(min(combined_theta_data$nucleotide_diversity)), 
                              #ceiling(max(combined_theta_data$nucleotide_diversity)), 
                              #by = 10)) +  # Set y-axis breaks every 10 units
  labs(title = "Nucleotide Diversity (π) across Populations", x = "Population", y = expression("log[10](θπ)")) +
  custom_theme

# Save the plots
ggsave("results/plots/hap2/theta/tajimasD_violinplot.tiff", plot = p_tajimasD_violin, width = 6, height = 6)
ggsave("results/plots/hap2/theta/nucleotide_diversity_violinplot.tiff", plot = p_nucleotide_div_violin, width = 6, height = 6)

# Split Violin Plot with Boxplot Overlay using group colors for Core vs. Edge.
# Define a two-color mapping for groups based on one representative color for each group.
#group_colors <- c("Core" = population_colors["C.OH.1"], "Edge" = population_colors["E.ON.5"])

#split_violin_plot <- ggplot(non_zero_data, aes(x = population, y = tajimasD/sites, fill = population_group)) +
  #geom_violin(trim = FALSE, position = position_dodge(0.9)) +
  #geom_boxplot(width = 0.2, outlier.shape = NA, position = position_dodge(0.9), color = "black") +
  #scale_fill_manual(values = group_colors) +
  #labs(title = "Split Violin and Boxplot of Tajima's D (Excluding Zeros)", x = "Population", y = "Tajima's D") +
  #custom_theme
#ggsave("results/plots/hap2/theta/tajimasD_split_violin_boxplot.tiff", plot = split_violin_plot, width = 6, height = 6)