# Loading necessary libraries
library(ggplot2)

# Process LR files
args <- commandArgs(trailingOnly = TRUE)
print(args)
input_file <- args[1]
pop_names <- args[2]
plot_file <- args[3]


# Read population names and covariance matrix
pop <- read.table(pop_names, header = FALSE, sep = "\t")
C <- as.matrix(read.table(input_file, header = FALSE))
e <- eigen(C)


# Prepare data for ggplot
pca_data <- data.frame(PC1 = e$vectors[, 1], PC2 = e$vectors[, 2], Population = pop[, 1])


cluster_colors <- c(
  # Existing colors:
  "C.OH.1" = "#DA70D6",    # orchid
  "C.MI.2" = "#9932CC",    # darkorchid
  "C.ON.3" = "#551A8B",    # purple4
  "C.IN.2" = "#FFA07A",    # lightsalmon1
  "E.ON.5" = "#00FFFF",    # cyan
  "E.ON.8" = "#7FFFD4",    # aquamarine
  "E.ON.9" = "#00868B",    # turquoise4
  
  # New core colors (4):
  "C.IL.1" = "#9370DB",    # mediumpurple
  "C.MN.4" = "#BA55D3",    # mediumorchid
  "C.WI.1" = "#663399",    # rebeccapurple
  "C.WI.2" = "#9400D3",    # darkviolet
  
  # New edge colors (4):
  "E.MI.7" = "#20B2AA",    # lightseagreen
  "E.MN.2" = "#48D1CC",    # mediumturquoise
  "E.MN.3" = "#5F9EA0",    # cadetblue
  "E.NH.1" = "#008B8B",    # darkcyan
  "E.NY.1" = "#00CED1",    # dark turqoise 
  "E.WI.4" = "#40E0D0",    # turquoise
  "E.WI.9" = "#009999"     # balanced teal
)


# Set the Population column as a factor with specific levels to maintain consistent coloring
pca_data$Population <- factor(pca_data$Population, levels = names(cluster_colors))

# Define a custom theme
custom_theme <- theme(panel.background = element_rect(fill = "white", color = NA),
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black", linewidth = 1), 
                      axis.ticks = element_line(linewidth = 1), 
                      axis.ticks.length = unit(0.25, "cm"), 
                      plot.title = element_text(hjust = 0.5), 
                      axis.text = element_text(size=12),
                      axis.title = element_text(size=12),
                      legend.position = "top", 
                      legend.margin = margin(0, 0, 0, 0), 
                      legend.text = element_text(size=12), 
                      legend.title = element_text(size=10))

# Create the plot using ggplot
gg_pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Population)) +
  geom_point(size = 3) +
  labs(title = "Individual Allele Frequency", x = "PC1", y = "PC2") +
  scale_color_manual(values = cluster_colors, name = "Population") +
  custom_theme

# Save the plot as a PNG file
ggsave(plot_file, plot = gg_pca_plot, width = 10, height = 8)
