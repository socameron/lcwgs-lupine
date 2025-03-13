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


# Initialize an empty data frame to combine data for all populations
combined_theta_data <- data.frame()


# Loop through populations and read their respective theta files
for (pop in populations) {
  file_path <- paste0("results/theta/hap2/", pop, "_out.thetasWindow.gz.pestPG")
  
  # Read the theta data
  theta_data <- read_theta_data(file_path)
  
  # Replace population names with new ones
  pop <- sapply(pop, function(x) {
    idx <- match(x, old_levels)
    if (!is.na(idx)) new_levels[idx] else x
  })
  
  # Add population identifier and group (Core/Edge)
  theta_data <- theta_data %>%
    mutate(population = pop,
           population_group = ifelse(pop %in% core_populations, "Core", "Edge"))
  
  # Combine data for all populations
  combined_theta_data <- rbind(combined_theta_data, theta_data)
}

# Reorder populations by core and edge
combined_theta_data$population <- factor(combined_theta_data$population, 
                                         levels = c("C.IN.2", "C.OH.1", "C.MI.2", "C.ON.1", "C.ON.3", 
                                                    "E.MI.3", "E.ON.5", "E.ON.8", "E.ON.9", "E.ON.12"))

# Count the number of Edge (E) and Core (C) populations
num_e_levels <- sum(grepl("^E", levels(combined_theta_data$population)))
num_c_levels <- sum(grepl("^C", levels(combined_theta_data$population)))

# Generate color palettes for Edge and Core populations
colours_e <- colorRampPalette(c("aquamarine", "cyan", "turquoise4", "dodgerblue"))(num_e_levels)
colours_c <- colorRampPalette(c("lightsalmon1", "orchid", "darkorchid", "purple4"))(num_c_levels)

# Assign colors to each population
population_colors <- vector("character", length(levels(combined_theta_data$population)))
names(population_colors) <- levels(combined_theta_data$population)


for (level in levels(combined_theta_data$population)) {
  if (startsWith(level, "E")) {
    population_colors[level] <- colours_e[which(grepl(level, levels(combined_theta_data$population)[grepl("^E", levels(combined_theta_data$population))]))]
  } else if (startsWith(level, "C")) {
    population_colors[level] <- colours_c[which(grepl(level, levels(combined_theta_data$population)[grepl("^C", levels(combined_theta_data$population))]))]
  }
}

# Create a color mapping for core and edge populations
#population_colors <- c("darkorchid3", "darkorchid3", "darkorchid3", "darkorchid3", # Core colors
                       #"cyan3", "cyan3", "cyan3") # Edge colors

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

mean_values <- combined_theta_data %>%
  group_by(population) %>%
  summarise(mean_diversity = mean(nucleotide_diversity/sites, na.rm = TRUE))

median_values <- combined_theta_data %>%
  group_by(population) %>%
  summarise(median_diversity = median(nucleotide_diversity/sites, na.rm = TRUE))


# Plotting boxplot for both Tajima's D and Nucleotide Diversity
p_tajimasD <- combined_theta_data %>%
  filter(!sites == 0) %>%
  ggplot(aes(x = population, y = tajimasD/sites, fill = population)) +
  geom_boxplot() +
  scale_fill_manual(values = population_colors) +
  labs(title = "Tajima's D across Populations", x = "Population", y = "Tajima's D") +
  custom_theme

p_nucleotide_div <- combined_theta_data %>%
  filter(!sites == 0) %>% 
  ggplot(aes(x = population, y = nucleotide_diversity/sites, fill = population)) +
  geom_boxplot() +
  geom_hline(data = mean_values, aes(yintercept = mean_diversity), color = "blue", linetype = "dashed", size = 0.5) +
  geom_hline(data = median_values, aes(yintercept = median_diversity), color = "red", linetype = "solid", size = 0.5) +
  scale_fill_manual(values = population_colors) +
  labs(title = "Nucleotide Diversity (π) across Populations", x = "Population", y = "Nucleotide Diversity (θπ)") +
  custom_theme

# Save the plots
ggsave("results/plots/hap2/theta/tajimasD_boxplot.tiff", plot = p_tajimasD, width = 10, height = 6)
ggsave("results/plots/hap2/theta/nucleotide_diversity_boxplot.tiff", plot = p_nucleotide_div, width = 10, height = 6)

# Calculate the 90th percentile of the nucleotide diversity
threshold <- quantile(combined_theta_data$nucleotide_diversity, 0.90, na.rm = TRUE)

# Plotting violin plot for both Tajima's D and Nucleotide Diversity
p_tajimasD_violin <- combined_theta_data %>%
  filter(!sites == 0) %>%
  ggplot(aes(x = population, y = tajimasD/sites, fill = population)) +
  geom_violin() + 
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") +
  scale_fill_manual(values = population_colors) +
  labs(title = "Tajima's D across Populations", x = "Population", y = "Tajima's D") +
  custom_theme

p_nucleotide_div_violin <- combined_theta_data %>%
  filter(!sites == 0) %>%
  ggplot(aes(x = population, y = log10(nucleotide_diversity/sites), fill = population)) +
  geom_violin(width=1) + 
  geom_hline(data = mean_values, aes(yintercept = mean_diversity), color = "blue", linetype = "dashed", size = 0.5) +
  geom_hline(data = median_values, aes(yintercept = median_diversity), color = "red", linetype = "solid", size = 0.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") +
  scale_fill_manual(values = population_colors) +
  #scale_y_continuous(breaks = seq(floor(min(combined_theta_data$nucleotide_diversity)), 
                              #ceiling(max(combined_theta_data$nucleotide_diversity)), 
                              #by = 10)) +  # Set y-axis breaks every 10 units
  labs(title = "Nucleotide Diversity (π) across Populations", x = "Population", y = expression("log[10](θπ)")) +
  custom_theme

# Save the plots
ggsave("results/plots/hap2/theta/tajimasD_violinplot.tiff", plot = p_tajimasD_violin, width = 6, height = 6)
ggsave("results/plots/hap2/theta/nucleotide_diversity_violinplot.tiff", plot = p_nucleotide_div_violin, width = 6, height = 6)
