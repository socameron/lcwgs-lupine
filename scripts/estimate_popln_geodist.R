# calc_popln_geodist.R
library(geosphere)
args <- commandArgs(trailingOnly = TRUE)

# Read the data from the CSV file given as the first argument
data <- read.csv(args[1], stringsAsFactors = FALSE)

# Function to calculate pairwise distances
calculate_distances <- function(df) {
  coords <- as.matrix(df[, c('Latitude', 'Longitude')])
  distances <- distm(coords, fun = distHaversine)
  distances_km <- distances / 1000
  return(distances_km)
}

# Calculate the distances
distance_matrix <- calculate_distances(data)

# Convert the distance matrix to a data frame
distance_df <- as.data.frame(as.table(distance_matrix))
names(distance_df) <- c('Sample1', 'Sample2', 'Distance_km')

# Write the distance matrix to a CSV file given as the second argument
write.csv(distance_df, args[2], row.names = FALSE)
