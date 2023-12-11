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

# Function to normalize
norm <- function(x) x/sum(x)

# Process the input file
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
plot_file <- args[2]

# Read data
sfs <- read_winSFS(input_file)
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
  theme_bw() +  # White background theme
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_blank(),      # Remove panel border
        axis.line = element_line(colour = "black")) +  # Add axis lines
  geom_text(aes(label = round(Proportions, 3)), vjust = -0.5, size = 3)

# Save the plot
ggsave(plot_file, plot = gg_sfs_plot, width = 10, height = 8)
