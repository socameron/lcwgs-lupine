# Load necessary libraries
library(ggplot2)
library(dplyr)

# Define the old population IDs (from your Snakemake wildcard)
old_populations <- c("HPW", "IDNP-MW", "LCTGP", "MFNP", "PPP", "RLPLV", "SWCP", "APB", "RSB", "CPB",
                      "BSL", "BSNA", "FMB", "GRAY", "NBWA", "NGP", "PBBT", "UWA")
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
combined_inbreeding_data <- data.frame()

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

# Loop through populations and read their respective inbreeding files
for (pop_old in populations) {
  file_path <- paste0("results/ngsF/hap2/by_popln/", pop_old, "_ngsF-HMM_inbreeding.indF")
  
  # Read the inbreeding data
  inbreeding_data <- read_inbreeding_data(file_path)$inbreeding_data

  # Convert old popln names to new names
  new_pop <- pop_mapping[pop_old]
  
  # Replace population names
  pop_group <- ifelse(startsWith(new_pop, "C."), "Core", "Edge")
  
  # Extract only the 'rab' column and add a population identifier
  inbreeding_data <- inbreeding_data %>%
    select(InbreedingCoefficient) %>%
    mutate(population = new_pop,
           population_group = pop_group)
  
  # Combine data for all populations
  combined_inbreeding_data <- rbind(combined_inbreeding_data, inbreeding_data)
}

# Reorder populations by core and edge
combined_inbreeding_data$population <- factor(combined_inbreeding_data$population, levels = c(
  "C.IN.2", "C.OH.1", "C.MI.2", "C.ON.3", "C.IL.1", "C.MN.4", "C.WI.1", "C.WI.2",
  "E.ON.5", "E.ON.8", "E.ON.9", "E.MI.7", "E.MN.2", "E.MN.3", "E.NH.1", "E.NY.1",  "E.WI.4", "E.WI.9"
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

# Create a boxplot for inbreeding coefficients
p_inbreeding <- combined_inbreeding_data %>%
  ggplot(aes(x = population, y = InbreedingCoefficient, fill = population)) +
  geom_boxplot() +
  scale_fill_manual(values = population_colors) +
  scale_x_discrete(labels = function(x) sub("^[CE]\\.", "", x)) + # remove C. and E. precursors of population names
  labs(title = "Inbreeding Coefficients across Populations", x = "Population", y = "Inbreeding Coefficient (F)") +
  custom_theme

# Save the boxplot
ggsave("results/plots/hap2/ngsF/ngsF-HMM_inbreeding_coeff_boxplot.tiff", plot = p_inbreeding, width = 8, height = 6)


inbreeding_summary <- combined_inbreeding_data %>%
  group_by(population) %>%
  summarise(
    mean_F = mean(InbreedingCoefficient, na.rm = TRUE),
    sd_F = sd(InbreedingCoefficient, na.rm = TRUE),
    n = n(),
    se_F = sd_F / sqrt(n)
  )

p_mean_inbreeding <- ggplot(inbreeding_summary, aes(x = population, y = mean_F, color = population)) +
  geom_point(size = 6) +
  geom_errorbar(aes(ymin = mean_F - se_F, ymax = mean_F + se_F), 
                width = 0.75, linewidth = 1.5) +
  scale_color_manual(values = population_colors) +
  scale_x_discrete(labels = function(x) sub("^[CE]\\.", "", x)) + # remove C. and E. precursors of population names
  labs(
    title = "Mean Inbreeding Coefficient by Population (F[IS])",
    x = "Population",
    y = expression("Inbreeding Coefficient (F)")
  ) +
  custom_theme

# Save the point plot
ggsave("results/plots/hap2/ngsF/ngsF-HMM_inbreeding_coeff_mean.tiff", plot = p_mean_inbreeding, width = 8, height = 6)
