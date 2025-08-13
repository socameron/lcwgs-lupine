#!/usr/bin/env Rscript
#
# scripts/plot_ngsAdmix_mapmixture.R
#

# 0) parse command‐line args
args        <- commandArgs(trailingOnly=TRUE)
qfile       <- args[1]   # e.g. results/ngsAdmix/.../canonical_K18.qopt
meta_csv    <- args[2]   # e.g. data/lists/.../all_samples_to_popln_bam_order.csv
coords_csv  <- args[3]   # e.g. data/lists/.../all_popln_geo_coord.csv
barplot_out <- args[4]   # e.g. results/plots/.../K18_admixture_barplot.png
piemap_out  <- args[5]   # e.g. results/plots/.../K18_pie_charts_map_mapmixture.png
combined_out<- args[6]   # e.g. results/plots/.../K18_combined_mapmixture.png

# Use when troubleshooting in R!
args        <- commandArgs(trailingOnly=TRUE)
qfile       <- "results/ngsAdmix/hap2/canonical_K18.qopt"
meta_csv    <- "data/lists/hap2/all_samples_to_popln_bam_order.csv"
coords_csv  <- "data/lists/hap2/all_popln_geo_coord.csv"
barplot_out <- "results/plots/hap2/ngsAdmix/K18_admixture_barplot_2.png"
piemap_out  <- "results/plots/hap2/ngsAdmix/K18_pie_charts_map_mapmixture.png"
combined_out<- "results/plots/hap2/ngsAdmix/K18_combined_mapmixture.png"


# 1) libraries
library(tidyverse)
library(RColorBrewer)
library(mapmixture)
library(ggpubr)
library(rworldmap)
library(rworldxtra)
library(rnaturalearthhires)
library(sf)
library(maps)

# 2) custom chart themes plus other map stuff
custom_theme <- theme_minimal() +
  theme(
    panel.grid     = element_blank(),
    axis.line      = element_line(colour="black"),
    axis.ticks     = element_line(),
    plot.title     = element_text(hjust=0.5),
    legend.position= "top"
  )


# 3) load inputs
Qmat   <- read.table(qfile,    header = FALSE, stringsAsFactors = FALSE)
meta   <- read.csv(meta_csv,   stringsAsFactors = FALSE)
coords <- read.csv(coords_csv, stringsAsFactors = FALSE)

# 4) build `admixture_df`: first col `site`, second `individual`, then cluster1…clusterK
K <- ncol(Qmat)
colnames(Qmat) <- paste0("cluster", seq_len(K))
admixture_df <- meta %>%
  transmute(
    site       = Site_code2,       # your population code
    individual = Sample_name        # your sample ID
  ) %>%
  bind_cols(Qmat)

# 5) build `coords_df`: must have site, lon, lat
coords_df <- coords %>%
  transmute(
    site = Site_code2,
    lat  = Latitude,
    lon  = Longitude
  )

# 6) download maps features
# note: use of Natural Earth Hires might require custom installation on R on the cluster, prior to automated use
world_all <- ne_countries(scale="large", returnclass="sf")[, "geometry"] %>% st_transform(crs_proj) #from world land 
lakes_all  <- ne_download(scale = "large", type = "lakes", category = "physical", returnclass="sf") #from Natural Earth
rivers_all <- ne_download(scale = "large", type = "rivers_lake_centerlines", category="physical", returnclass="sf") #from Natural Earth
oceans_all <- ne_download(scale = "large", type = "ocean", category = "physical", returnclass = "sf")
states_all <- ne_states(country = c("United States of America", "Canada"), returnclass = "sf")
lands_all <- ne_countries(scale="medium", returnclass="sf")

# 7) create a colour palette for the K clusters
cluster_pal   <- colorRampPalette(brewer.pal(8, "Set2"))(K)
cluster_names <- paste0("Cluster ", seq_len(K))

# 8) STRUCTURE‐style barplot (fill by Cluster, faceted by site)
qlong <- admixture_df %>%
  pivot_longer(
    cols      = paste0("cluster", seq_len(K)),
    names_to  = "Cluster",
    values_to = "Prop"
  )
# using ggplot
p_bar <- ggplot(qlong, aes(x = individual, y = Prop, fill = Cluster)) +
  geom_col(width = 1) +
  facet_grid(~ site, scales="free_x", space="free_x") +
  scale_fill_manual(values = cluster_pal, labels = cluster_names) +
  labs(
    title = paste0("Individual admixture (K=", K, ")"),
    x     = NULL,
    y     = "Proportion",
    fill  = "Cluster"
  ) +
  custom_theme +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(face = "bold")
  )

ggsave(barplot_out, p_bar, width = 10, height = 2.5, dpi = 300)

# using mapmixture built-in functions
structure_barplot <- structure_plot(
  admixture_df = admixture_df,
  type          = "structure",
  cluster_cols  = cluster_pal,
  cluster_names = cluster_names,   # optional, will be used in the legend
  site_dividers = TRUE,            # draw little vertical dividers between sites
  divider_width = 0.4,             # thickness of those dividers
  site_order    = unique(admixture_df$site),   # in the order your sites appear
  labels        = "site",          # show each site code under its block
  flip_axis     = FALSE,           # if TRUE, sites go on the y‐axis
  site_ticks_size = -0.05,         # how far labels sit from the axis
  site_labels_y   = -0.35,         # vertical position of the site labels
  site_labels_size= 2.2            # font size of those labels
) +
  custom_theme +
  theme(
    axis.title.y = element_text(size = 8, hjust = 1),
    axis.text.y  = element_text(size = 5)
  )

ggsave("results/plots/hap2/ngsAdmix/K18_structure_plot.png", structure_barplot, width = 10, height = 2.5, dpi = 300)



# 9) Admixture Graphs - using mapmixture pie graphs


#### MAP 1 - Use 'ggplot2 functions ####
# We build a base map with whichever layers we use, then add the piecharts.
# This is because mapmixture only allows one sf layer as the base map and it's difficult to add additional layers
admix_coords <- meta %>%
  transmute(site       = Site_code2,
            individual = Sample_name) %>%
  bind_cols(Qmat) %>% 
  left_join(coords %>% select(site = Site_code2, lat = Latitude, lon = Longitude), by = "site") %>%
  select(site, lat, lon, paste0("cluster", 1:K))
cluster_cols <- 4:ncol(admix_coords)

avg_admix_coords <- admix_coords %>%
  group_by(site, lat, lon) %>%
  summarize(across(starts_with("cluster"), mean), .groups="drop")

# convert coordinates to crs_proj (5070)
crs_proj <- 5070 # NAD 83 coordinate system
admix_sf <- st_as_sf(
  avg_admix_coords,
  coords = c("lon","lat"),
  crs    = 4326,
  remove = FALSE)
admix_sf <- st_transform(admix_sf, crs_proj)
coords_proj <- st_coordinates(admix_sf)
avg_admix_coords$lon <- coords_proj[,1]
avg_admix_coords$lat <- coords_proj[,2]

# convert bounding box to crs_proj
bb_ll    <- c(xmin=-100, xmax=-70, ymin=35, ymax=55)
bb_proj  <- transform_bbox(bb_ll, crs_proj)

world <- st_transform(world_all, crs_proj)
lakes  <- st_transform(lakes_all,  crs_proj)
rivers <- st_transform(rivers_all, crs_proj)
oceans <- st_transform(oceans_all, crs_proj)
states <- st_transform(states, crs_proj)

p_map <- ggplot() +
  geom_sf(data = oceans, fill = "#cbe8f6", color = NA) +
  geom_sf(data = world,   fill="#f0f0f0", color="grey50", size = 0.2) +
  geom_sf(data = lakes,   fill="#a6cee3", color=NA) +
  #geom_sf(data = rivers,  color="#BEE6FE", size=0.3) +BEE6FE
  geom_sf(data = states, fill = NA, color = "grey70", size = 0.2) +
  add_pie_charts(
    df           = avg_admix_coords,
    admix_columns= cluster_cols, 
    lat_column   = "lat",
    lon_column   = "lon",
    pie_colours  = cluster_pal,
    pie_size     = 0.7,
    border       = 0.2,
    opacity      = 1
  ) +
  coord_sf(
    xlim        = c(bb_proj["xmin"], bb_proj["xmax"]),
    ylim        = c(bb_proj["ymin"], bb_proj["ymax"]),
    expand      = FALSE
  ) +
  theme_void() +
  theme(legend.position="none") +
  labs(title = paste0("Mean admixture per site (K=",K,")"))

ggsave(piemap_out, p_map, width = 6, height = 6, dpi = 300)

# 10) combine into one figure
p_comb <- ggarrange(
  p_bar + theme(legend.position="none"),
  p_map,
  ncol   = 1,
  heights= c(0.25, 0.75)
)
ggsave(combined_out, p_comb, width = 10, height = 8, dpi = 300)







































#### ARCHIVE ####

#### MAP 1 - Use 'mapmixture' function ####
# NOTE: it is difficult to add additional layers after mapmixture() 

# create bounding box
crs_proj    <- 5070
bb_ll       <- c(xmin=-100, xmax=-70, ymin=35, ymax=55)
bb_proj     <- transform_bbox(bb_ll, crs_proj)

# create pie sizes based on coordinates (OPTIONAL)
xy <- coords_df %>% select(lon, lat)
dmat <- as.matrix(dist(xy))
dmin <- min(dmat[dmat>0])
# we shrink the pies to 30–40% of that
pie_size <- dmin * 0.35

# create base map that includes oceans, land, lakes, rivers, and state/provinvical boundaries
map_xlim <- c(-100, -70)
map_ylim <- c(35,  55)
base_map <- ggplot() +
  geom_sf(data = oceans_all, fill = "#cbe8f6", color = NA) +
  geom_sf(data = world_all,  fill = "#f0f0f0", color = NA) +
  geom_sf(data = lands_all,  fill = NA,        color = "grey50", size = 0.2) +
  geom_sf(data = lakes_all,  fill = "#a6cee3", color = NA) +
  geom_sf(data = rivers_all, color = "#9ecae1", size = 0.3) +
  geom_sf(data = states_all, fill = NA,        color = "grey70", size = 0.2) +
  coord_sf(
    crs    = crs_proj,
    xlim   = map_xlim,
    ylim   = map_ylim,
    expand = FALSE
  ) +
  theme_void()



# use mapmixture to create graph
p_map <- mapmixture(
  admixture_df  = admixture_df,
  basemap        = world_all,
  coords_df      = coords_df,
  cluster_cols   = cluster_pal,
  cluster_names  = cluster_names,
  crs            = crs_proj,
  pie_size       = 1,
  pie_border     = 0.2,
  pie_border_col = "black",
  pie_opacity    = 1,
  land_colour    = "#f0f0f0",
  sea_colour     = "#cbe8f6",
  arrow          = FALSE,
  scalebar       = FALSE,
  boundary       = c(xmin=-100, xmax=-70, ymin=35, ymax=55)) +
  geom_sf(data = lakes_all,  fill = "#a6cee3", color = NA, inherit.aes = FALSE) +
  geom_sf(data = rivers_all, color = "#9ecae1", size = 0.3, inherit.aes = FALSE) +
  geom_sf(data = states_all, fill = NA,        color = "grey70", size = 0.2, inherit.aes = FALSE) +
  theme_void() +
  theme(legend.position="none",
        plot.title      = element_text(hjust=0.5)) +
  labs(title = paste0("Mean admixture per site (K=", K, ")"))

ggsave(piemap_out, p_map, width = 6, height = 6, dpi = 300)
