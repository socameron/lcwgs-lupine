#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(ggrepel)

args      <- commandArgs(trailingOnly = TRUE)
cov_file  <- args[1]
meta_csv  <- args[2]
out_png   <- args[3]

# 1) Read metadata and covariance
meta <- read.csv(meta_csv, stringsAsFactors = FALSE)
C    <- as.matrix(read.table(cov_file, header = FALSE))
e    <- eigen(C)

# 2) Build PCA data.frame
pca_data <- data.frame(
  PC1        = e$vectors[,1],
  PC2        = e$vectors[,2],
  Population = meta$Site_code2
)

# 3) Color palette & factor
cluster_colors <- c(
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
pca_data$Population <- factor(pca_data$Population, levels=names(cluster_colors))

# 4) Compute centroids
centroids <- pca_data %>%
  group_by(Population) %>%
  summarize(x=mean(PC1), y=mean(PC2), .groups="drop")

# 5) Plot with ggrepel
custom_theme <- theme(panel.border = element_blank(), 
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black", linewidth = 1), 
                      axis.ticks = element_line(linewidth = 1), 
                      axis.ticks.length = unit(0.25, "cm"), 
                      plot.title = element_text(hjust = 0.5), 
                      axis.text = element_text(size=12, angle = 45, hjust = 1), 
                      axis.title = element_text(size=12),
                      legend.position = "none", 
                      legend.margin = margin(0, 0, 0, 0), 
                      legend.text = element_text(size=12), 
                      legend.title = element_text(size=10))

ggplot(pca_data, aes(x=PC1, y=PC2, color=Population)) +
  geom_point(size=3, alpha=0.7) +
  scale_color_manual(values=cluster_colors) +
  scale_x_reverse() + # flip x-axis so 'western' populations are on the left. 
  labs(title="Individual PCA", x="PC1", y="PC2") +
  geom_text_repel(
    data          = centroids,
    aes(x=x, y=y, label=Population, color=Population),
    fontface      = "bold",
    size          = 5,
    show.legend   = FALSE,
    nudge_y       = 0.05,
    box.padding   = 0.05,
    point.padding = 0.05,
    segment.color = "grey50",
    segment.size  = 0.2) +
  custom_theme

ggsave(out_png, width=10, height=8)

# Archived geom_text with lines


#geom_text_repel(
  #data          = centroids,
  #aes(x = x, y = y, label = Population, color = Population),
  #fontface      = "bold",
  #size          = 5,
  #show.legend   = FALSE,
  #direction     = "both",                    # allow movement in x & y
  #force         = 2,                         # stronger repulsion
  #box.padding   = unit(0.6, "lines"),        # gap around the text box
  #point.padding = unit(0.6, "lines"),        # gap around the data point
  #segment.color = "grey50",                  # little leader line
  #segment.size  = 0.5,
  #max.overlaps  = Inf                        # never drop a label
  #)