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

# Function to read inbreeding data from .indF files
read_inbreeding_data <- function(file_path) {
  # Read the .indF file, skipping the first line (log-likelihood)
  all_data <- read.table(file_path, header = FALSE, sep = "\t", fill = TRUE, skip=1)
  
  # Dynamically identify where the individual inbreeding coefficients end
  # Assuming the MAF section has only one column, we find the break
  maf_start <- which(is.na(all_data$V2))[1]
  
  # Extract inbreeding coefficients and transition rates
  inbreeding_data <- all_data[1:(maf_start - 1), 1:2]
  
  # Rename columns
  colnames(inbreeding_data) <- c("InbreedingCoefficient", "TransitionRate")
  
  # Minor allele frequencies (after individual inbreeding coefficients)
  maf_data <- all_data[maf_start:nrow(all_data), 1]
  
  return(list(inbreeding_data = inbreeding_data, maf_data = maf_data))
}

# Initialize an empty data frame to combine data for all populations
combined_inbreeding_data <- data.frame()

# Loop through populations and read their respective inbreeding files
for (pop in populations) {
  file_path <- paste0("results/ngsF/hap2/by_popln/", pop, "_ngsF-HMM_inbreeding.indF")
  
  # Read the inbreeding data
  inbreeding_data <- read_inbreeding_data(file_path)$inbreeding_data  # Only get inbreeding data
  
  # Replace population names with new ones
  pop <- new_levels[match(pop, old_levels)]
  
  # Add population identifier and group (Core/Edge)
  inbreeding_data <- inbreeding_data %>%
    mutate(population = pop,
           population_group = ifelse(pop %in% core_populations, "Core", "Edge"))
  
  # Combine data for all populations
  combined_inbreeding_data <- rbind(combined_inbreeding_data, inbreeding_data)
}

# Reorder populations by core and edge
combined_inbreeding_data$population <- factor(combined_inbreeding_data$population, 
                                              levels = c("C.IN.2", "C.OH.1", "C.MI.2", "C.ON.1", "C.ON.3", 
                                                         "E.MI.3", "E.ON.5", "E.ON.8", "E.ON.9", "E.ON.12"))

# Count the number of Edge (E) and Core (C) populations
num_e_levels <- sum(grepl("^E", levels(combined_inbreeding_data$population)))
num_c_levels <- sum(grepl("^C", levels(combined_inbreeding_data$population)))

# Generate color palettes for Edge and Core populations
colours_e <- colorRampPalette(c("aquamarine", "cyan", "turquoise4", "dodgerblue"))(num_e_levels)
colours_c <- colorRampPalette(c("lightsalmon1", "orchid", "darkorchid", "purple4"))(num_c_levels)

# Assign colors to each population
population_colors <- setNames(c(colours_c, colours_e), levels(combined_inbreeding_data$population))

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

# Create a boxplot for inbreeding coefficients
p_inbreeding <- combined_inbreeding_data %>%
  ggplot(aes(x = population, y = InbreedingCoefficient, fill = population)) +
  geom_boxplot() +
  scale_fill_manual(values = population_colors) +
  labs(title = "Inbreeding Coefficients across Populations", x = "Population", y = "Inbreeding Coefficient") +
  custom_theme

# Save the boxplot
ggsave("results/plots/hap2/ngsF/ngsF-HMM_inbreeding_coeff.tiff", plot = p_inbreeding, width = 6, height = 6)
