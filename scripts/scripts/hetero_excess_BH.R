args <- commandArgs(trailingOnly = TRUE)
input_files <- args[-length(args)]  # All but last argument
plot_file <- tail(args, n=1)  # Last argument

library(ggplot2)
library(gridExtra)

plot_list <- list()

# Loop through each file and create a histogram
for (file in input_files) {
  tryCatch({
    # Extract population and BH_VAR from file name
    file_parts <- unlist(strsplit(basename(file), "_"))
    population <- file_parts[2]
    BH_VAR <- sub("BH", "", sub(".lr", "", file_parts[length(file_parts)]))

    hwe_data <- read.table(gzfile(file), header = FALSE, sep="\t")
    p <- ggplot(hwe_data, aes(x=V7)) + 
      geom_histogram(binwidth = 0.05, fill="blue", color="black") +
      labs(title=paste(population, ": Benjamini-Hochberg crit =", BH_VAR), x="F value", y="Frequency")
    plot_list[[length(plot_list) + 1]] <- p
  }, error = function(e) {
    cat("Error in file:", file, "\n", e$message, "\n")
  })
}

# Combine the plots
combined_plot <- do.call(grid.arrange, c(plot_list, ncol=3))

# Save the combined plot
ggsave(filename = plot_file, plot = combined_plot, device = "png", width = 48, height = 25)
