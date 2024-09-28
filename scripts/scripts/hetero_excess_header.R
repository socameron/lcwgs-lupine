args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
plot_file <- args[2]

# Read the gzipped HWE file
hwe_data <- read.table(gzfile(input_file), header = TRUE, sep="\t")

# Load ggplot2
library(ggplot2)
hwe_data$F <- as.numeric(hwe_data$F)

# Create the histogram
p <- ggplot(hwe_data, aes(x=F)) + 
  geom_histogram(binwidth = 0.05, fill="red", color="black") +
  labs(title="Histogram of F Values", x="F value", y="Frequency")

# Save the plot
ggsave(plot_file, plot = p, width = 10, height = 8)
