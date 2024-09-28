#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(dplyr)
library(minpack.lm)  # For nonlinear least squares fitting

# Function to plot LD for a given population
plot_ld <- function(ld_file, output_file) {
  # Read the .ld file
  ld_data <- read.table(ld_file, header=TRUE, stringsAsFactors=FALSE)
  
  # Calculate distance
  ld_data <- ld_data %>%
    mutate(Distance = abs(BP_A - BP_B)) %>%
    sample_n(1000000)
  
  # Fit the Hill and Robertson model
  hill_robertson_model <- function(d, N) {
    return(1 / (1 + 4 * N * 1e-8 * d))
  }
  
  # Initial guess for N
  initial_guess <- list(N = 1e4)
  
  # Fit the model
  fit <- nlsLM(R2 ~ hill_robertson_model(Distance, N), data = ld_data, start = initial_guess)
  
  # Get the fitted values
  fitted_values <- predict(fit)
  
  # Create a ggplot object with loess smoothing and theoretical LD decay curve
  p <- ggplot(ld_data, aes(x=Distance, y=R2)) +
    geom_point(alpha=0.3) +
    geom_smooth(method="loess", se=FALSE, color="blue", size=1) +
    geom_line(aes(y=fitted_values), color="red", size=1) +
    theme_minimal() +
    labs(title="Linkage Disequilibrium (LD)",
         x="Distance (bp)",
         y="R²") +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Save the plot as a JPEG
  ggsave(output_file, plot=p, width=10, height=8, dpi=300, device="jpeg")
}

# Main function to iterate over populations
main <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  ld_file <- args[1]
  output_file <- args[2]
  plot_ld(ld_file, output_file)
}

main()
