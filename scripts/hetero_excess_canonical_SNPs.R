args <- commandArgs(trailingOnly = TRUE)

# Read the gzipped HWE file
hwe_data <- read.table(gzfile(args[1]), header = FALSE, sep="\t")

# Load ggplot2
library(ggplot2)

# Create the histogram
p <- ggplot(hwe_data, aes(x=$V6)) + 
  geom_histogram(binwidth = 0.05, fill="blue", color="black") +
  labs(title="Histogram of F Values", x="F value", y="Frequency")

# Save the plot
ggsave(args[2], plot = p, width = 10, height = 8)
