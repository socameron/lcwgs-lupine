# Load necessary libraries
library(dplyr)

# Read command line arguments from Snakemake
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Read the ngsRelate output file
ngsrelate_data <- read.table(input_file, header = TRUE)

# Define a threshold for identifying clones (e.g., rab > 0.99)
clone_threshold <- 0.99

# Filter for pairs where relatedness (rab) is above the threshold
potential_clones <- ngsrelate_data %>%
  filter(rab > clone_threshold)

# Write the clone summary to a text file
if(nrow(potential_clones) > 0) {
  write.table(potential_clones, file = output_file, row.names = FALSE, col.names = TRUE, quote = FALSE)
} else {
  write("No clones detected", file = output_file)
}
