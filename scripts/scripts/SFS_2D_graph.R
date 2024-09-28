# Load necessary libraries
library(ggplot2)

# Function to read winSFS data
read_winSFS <- function(file) {
    lines <- readLines(file)
    shape_info <- strsplit(substring(lines[1], 8), "/")[[1]]
    nrow <- as.integer(shape_info[1])
    ncol <- as.integer(shape_info[2])
    data <- as.numeric(unlist(strsplit(lines[2], " ")))
    matrix(data, nrow = nrow, ncol = ncol, byrow = TRUE)
}

# Function to normalize
norm <- function(x) x/sum(x)

# Process the input file
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
plot_file <- args[2]

# Read and reshape data
sfs_matrix <- read_winSFS(input_file)

# Here you can add processing steps for the matrix
# For example, sum across rows or columns as needed

# Example: Sum across rows and then normalize
sfs <- apply(sfs_matrix, 1, sum)
sfs <- norm(sfs)

# Create a data frame for ggplot
sfs_data <- data.frame(Index = 1:length(sfs), Proportions = sfs)

# Extract filename for plot title
filename <- basename(input_file)

# Create the plot
gg_sfs_plot <- ggplot(sfs_data, aes(x = Index, y = Proportions)) +
  geom_bar(stat = "identity", fill = 'blue') +
  ggtitle(paste("winSFS Plot -", filename)) +
  xlab("Index") +
  ylab("Proportions") +
  theme_minimal() +
  geom_text(aes(label = round(Proportions, 3)), vjust = -0.5, size = 3)

# Save the plot
ggsave(plot_file, plot = gg_sfs_plot, width = 10, height = 8)
