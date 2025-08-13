#!/usr/bin/env Rscript

# Load arguments
args          <- commandArgs(trailingOnly = TRUE)
qfile         <- args[1]
metadata_file <- args[2]
coords_file   <- args[3]
barplot_out   <- args[4]
pie_out       <- args[5]
combined_out  <- args[6]

# Libraries
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(sf)
library(rworldmap)
library(rworldxtra)
library(ggspatial)     # for annotation_north_arrow + annotation_scale
library(grid)
library(mapmixture)

# — your custom colors & theme —
population_colors <- c(
  # Core populations:
  "C.OH.1" = "#DA70D6",    "C.MI.2" = "#9932CC",
  "C.ON.3" = "#551A8B",    "C.IN.2" = "#B19CD9",
  "C.IL.1" = "#9370DB",    "C.MN.4" = "#BA55D3",
  "C.WI.1" = "#663399",    "C.WI.2" = "#9400D3",
  # Edge populations:
  "E.ON.5" = "#00FFFF",    "E.ON.8" = "#7FFFD4",
  "E.ON.9" = "#00868B",    "E.MI.7" = "#20B2AA",
  "E.MN.2" = "#48D1CC",    "E.MN.3" = "#5F9EA0",
  "E.NH.1" = "#008B8B",    "E.NY.1" = "#00CED1",
  "E.WI.4" = "#40E0D0",    "E.WI.9" = "#009999"
)
custom_theme <- theme(
  panel.border      = element_blank(),
  panel.grid.major  = element_blank(),
  panel.grid.minor  = element_blank(),
  axis.line         = element_line(colour="black", size=0.5),
  axis.ticks        = element_line(size=0.5),
  axis.ticks.length = unit(0.25, "cm"),
  plot.title        = element_text(hjust=0.5),
  axis.text         = element_text(size=14, angle=45, hjust=1),
  axis.title        = element_text(size=14),
  legend.position   = "none",
  legend.margin     = margin(0,0,0,0),
  legend.text       = element_text(size=14),
  legend.title      = element_text(size=10)
)

# Read ancestry proportions
qmatrix <- read.table(qfile)
colnames(qmatrix) <- paste0("Cluster", seq_len(ncol(qmatrix)))

# Load sample metadata
metadata <- read.csv(metadata_file, stringsAsFactors=FALSE)
qmatrix$Ind  <- metadata$Ind
qmatrix$Site <- metadata$Site_code2

# Melt for barplot
qlong <- melt(qmatrix, id.vars=c("Ind","Site"))

# Use your custom palette for the K clusters:
cols <- rep(unname(population_colors), length.out=length(unique(qlong$variable)))

# — Barplot —
admix_bar <- ggplot(qlong, aes(x=Ind, y=value, fill=variable)) +
  geom_col() +
  facet_wrap(~Site, scales="free", ncol=3) +
  scale_fill_manual(values=cols) +
  ylab("Admixture proportion") +
  custom_theme +
  theme(
    axis.text.x     = element_blank(),
    axis.title.x    = element_blank(),
    strip.text      = element_text(size=10),
    legend.position = "top"
  )
ggsave(barplot_out, admix_bar, width=10, height=6, dpi=300)

# — Pie charts —
clusters  <- grep("^Cluster", names(qmatrix))
avg_admix <- aggregate(qmatrix[,clusters],
                       by=list(Site=qmatrix$Site), mean)
avg_admix <- melt(avg_admix, id.vars="Site")

pie_charts <- function(df, site){
  ggplot(subset(df, Site==site),
         aes(x="", y=value, fill=variable)) +
    geom_col(color="black") +
    coord_polar(theta="y") +
    scale_fill_manual(values=cols) +
    custom_theme +
    theme(
      legend.position="none",
      plot.margin=unit(rep(0.1,4),"cm")
    )
}

sites <- unique(avg_admix$Site)
pies  <- lapply(sites, pie_charts, df=avg_admix)
names(pies) <- sites

# Load coordinates
coords_all <- read.csv(coords_file)
coords <- coords_all %>%
  select(Site=Site_code, Lon=Longitude, Lat=Latitude) %>%
  filter(Site %in% sites) %>%
  arrange(match(Site, sites))

# Base map
map.outline <- getMap(resolution="high") %>%
  crop(extent(-180,180,-90,90)) %>%
  fortify()
basemap <- ggplot() +
  geom_polygon(data=map.outline,
               aes(x=long,y=lat,group=group),
               fill="grey85",color="black") +
  coord_quickmap() +
  custom_theme +
  # replace ggsn::north + scalebar with ggspatial:
  annotation_north_arrow(
    location = "tl", which_north = "true",
    style = north_arrow_fancy_orienteering()
  ) +
  annotation_scale(
    location = "bl", width_hint = 0.2
  )

# Place pies
radius <- 1
pie_layers <- mapply(function(p,coord){
  annotation_custom(
    grob=ggplotGrob(p),
    xmin=coord$Lon-radius, xmax=coord$Lon+radius,
    ymin=coord$Lat-radius, ymax=coord$Lat+radius
  )
}, pies, split(coords, coords$Site), SIMPLIFY=FALSE)

map_with_pies <- basemap + pie_layers
ggsave(pie_out, map_with_pies, width=10, height=8, dpi=300)

# — Combined —
combined <- ggarrange(
  admix_bar + labs(title="Individual Admixture (K=18)"),
  map_with_pies + labs(title="Mean admixture per site"),
  ncol=1
)
ggsave(combined_out, combined, width=12, height=12, dpi=300)
