# Load necessary libraries
library(ggplot2)

# Function to read winSFS data (1D)
read_winSFS <- function(file) {
    lines <- readLines(file)
    data <- as.numeric(unlist(strsplit(lines[2], " ")))
    
    # Check for any NA values in the data
    if (any(is.na(data))) {
        stop("NA values found in data conversion.")
    }
    
    return(data)
}

# Function to normalize excluding the first element
norm_excluding_first <- function(x) {
    x[-1] / sum(x[-1])
}

# Process the input file
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
plot_file <- args[2]

# Read data
sfs <- read_winSFS(input_file)

# Normalize data excluding the first data point (assuming it's the 0 class)
sfs_normalized <- norm_excluding_first(sfs)

# Create a data frame for ggplot
sfs_data <- data.frame(Index = 2:length(sfs), Proportions = sfs_normalized)

# Extract filename for plot title
filename <- basename(input_file)

# Create the plot
gg_sfs_plot <- ggplot(sfs_data, aes(x = Index, y = Proportions)) +
  geom_bar(stat = "identity", fill = 'blue') +
  ggtitle(paste("winSFS Plot -", filename)) +
  xlab("Index") +
  ylab("Proportions") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black")) +
  geom_text(aes(label = round(Proportions, 3)), vjust = -0.5, size = 3)

# Save the plot
ggsave(plot_file, plot = gg_sfs_plot, width = 10, height = 8)
