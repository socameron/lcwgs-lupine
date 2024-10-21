# Load necessary libraries
library(dplyr)
library(ggplot2)
library(geosphere)

# 1. Load population coordinates data
pop_coords <- read.csv("data/lists/batch_1_popln_only_geo_coord.csv")

# 2. Function to calculate distances between populations based on lat/long
calculate_distances <- function(pop1, pop2, coords) {
  # Retrieve coordinates for the first population using Site_code
  coord1 <- coords %>% filter(Site_code == pop1) %>% select(Longitude, Latitude) %>% unlist()
  
  # Retrieve coordinates for the second population using Site_code
  coord2 <- coords %>% filter(Site_code == pop2) %>% select(Longitude, Latitude) %>% unlist()
  
  # Check if both coordinates are valid (should have a length of 2)
  if (length(coord1) != 2 | length(coord2) != 2) {
    stop(paste("Invalid coordinates for", pop1, "or", pop2))
  }
  
  # Calculate geographic distance using the haversine formula (in kilometers)
  dist <- distHaversine(coord1, coord2) / 1000  # Distance in km
  return(dist)
}

# 3. Load Fst estimates
fst_files <- list.files(path = "results/realSFS/hap2/fst/", pattern = "*_fst_global.txt", full.names = TRUE)

# Function to extract weighted Fst and population pairs from filenames
get_fst_data <- function(file) {
  fst_values <- read.table(file, header = FALSE)
  weighted_fst <- fst_values$V2  # Second column is weighted Fst
  
  # Extract population pair from filename
  pop_pair <- basename(file) %>% gsub("_fst_global.txt", "", .)
  pops <- strsplit(pop_pair, "_")[[1]]
  pop1 <- pops[1]
  pop2 <- pops[2]
  
  # Return data frame with population pairs and Fst
  return(data.frame(Pop1 = pop1, Pop2 = pop2, Fst = weighted_fst))
}

# Combine Fst data from all files
fst_data <- do.call(rbind, lapply(fst_files, get_fst_data))

# 4. Remove same-population pairs (e.g., PPP-PPP)
fst_data <- fst_data %>% filter(Pop1 != Pop2)

# 5. Map Pop1 and Pop2 to their Site_code in pop_coords to use for distance calculations
fst_data <- fst_data %>%
  left_join(pop_coords %>% select(Site_code, New_code), by = c("Pop1" = "Site_code")) %>%
  rename(Pop1_New_code = New_code) %>%
  left_join(pop_coords %>% select(Site_code, New_code), by = c("Pop2" = "Site_code")) %>%
  rename(Pop2_New_code = New_code)

# 6. Create a unique identifier for each population pair, sorted alphabetically
fst_data <- fst_data %>%
  rowwise() %>%
  mutate(PairID = paste(sort(c(Pop1, Pop2)), collapse = "_"))

# 7. Remove duplicates based on the unique PairID
fst_data <- fst_data %>% distinct(PairID, .keep_all = TRUE)

# 8. Calculate distances between populations using Site_code and add to the dataframe
fst_data <- fst_data %>%
  rowwise() %>%
  mutate(Distance = calculate_distances(Pop1, Pop2, pop_coords))

# 9. Determine population pair type (Core-Core, Core-Edge, or Edge-Edge)
fst_data <- fst_data %>%
  left_join(pop_coords %>% select(New_code, Position), by = c("Pop1_New_code" = "New_code")) %>%
  rename(Pop1_Position = Position) %>%
  left_join(pop_coords %>% select(New_code, Position), by = c("Pop2_New_code" = "New_code")) %>%
  rename(Pop2_Position = Position) %>%
  mutate(PairType = case_when(
    Pop1_Position == "Core" & Pop2_Position == "Core" ~ "Core-Core",
    Pop1_Position == "Edge" & Pop2_Position == "Edge" ~ "Edge-Edge",
    TRUE ~ "Edge-Core"
  ))

# 10. Plot Fst vs. Distance with color coding for population pair type
# Custom theme settings
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

ggplot(fst_data, aes(x = Distance, y = Fst, color = PairType)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 1) +
  scale_color_manual(values = c("Edge-Core" = "#E3AF41", "Core-Core" = "#9932CC", "Edge-Edge" = "00868B")) +
  labs(title = "Fst vs. Geographic Distance", x = "Geographic Distance (km)", y = "Weighted Fst") +
  custom_theme +
  theme(legend.title = element_blank(), legend.position = "top")


# Save the plot
ggsave("results/plots/hap2/fst/fst_vs_distance.png", width = 8, height = 6)
