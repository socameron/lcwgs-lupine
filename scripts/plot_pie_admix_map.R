#!/usr/bin/env Rscript

# Load arguments
args <- commandArgs(trailingOnly = TRUE)
qfile <- args[1]
metadata_file <- args[2]
coords_file <- args[3]
barplot_out <- args[4]
pie_out <- args[5]
combined_out <- args[6]

# Libraries
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(sf)
library(rworldmap)
library(rworldxtra)
library(ggsn)
library(grid)

# Read Q-matrix
qmatrix <- read.table(qfile)
colnames(qmatrix) <- paste0("Cluster", 1:ncol(qmatrix))

# Load sample metadata
metadata <- read.csv(metadata_file)
qmatrix$Ind <- metadata$Ind
qmatrix$Site <- metadata$Site

# Melt for barplot
qlong <- melt(qmatrix, id.vars = c("Ind", "Site"))

# Colour palette
cols <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(unique(qlong$variable)))

# Barplot
admix_bar <- ggplot(qlong, aes(x = Ind, y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Site, scales = "free", ncol = 3) +
  scale_fill_manual(values = cols) +
  ylab("Admixture proportion") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        strip.text = element_text(size = 10), legend.position = "top")
ggsave(barplot_out, plot = admix_bar, width = 10, height = 6)

# Average admixture per site
clusters <- grep("Cluster", names(qmatrix))
avg_admix <- aggregate(qmatrix[, clusters], list(Site = qmatrix$Site), mean)
avg_admix <- melt(avg_admix, id.vars = "Site")

# Create pies
pie_charts <- function(df, site, cols) {
  ggplot(data = subset(df, Site == site), aes(x = "", y = value, fill = variable)) +
    geom_bar(stat = "identity", width = 1, color = "black") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = cols) +
    theme_void()
}

sites <- unique(avg_admix$Site)
pies <- lapply(sites, function(site) pie_charts(avg_admix, site, cols))
names(pies) <- sites

# Load coordinate data
coords_all <- read.csv(coords_file)

# Subset and rename to match site labels used in admixture results
coords <- coords_all %>%
  select(Site = Site_code, Lon = Longitude, Lat = Latitude) %>%
  filter(Site %in% sites) %>%
  arrange(match(Site, sites))


# Map base
map.outline <- crop(getMap(resolution = "high"), extent(-180, 180, -90, 90)) %>% fortify()
basemap <- ggplot() +
  geom_polygon(data = map.outline, aes(x = long, y = lat, group = group), fill = "grey85", color = "black") +
  coord_quickmap() +
  theme_void()

# Add pies
radius <- 1
pie_layers <- mapply(function(p, coord) {
  annotation_custom(grob = ggplotGrob(p),
                    xmin = coord$Lon - radius,
                    xmax = coord$Lon + radius,
                    ymin = coord$Lat - radius,
                    ymax = coord$Lat + radius)
}, pies, split(coords, coords$Site), SIMPLIFY = FALSE)

map_with_pies <- basemap + pie_layers
ggsave(pie_out, map_with_pies, width = 10, height = 8)

# Combined
combined <- ggarrange(admix_bar + labs(title = "Admixture barplot", tag = "A"),
                      map_with_pies + labs(title = "Mean admixture per site", tag = "B"),
                      ncol = 1)
ggsave(combined_out, combined, width = 12, height = 12)
