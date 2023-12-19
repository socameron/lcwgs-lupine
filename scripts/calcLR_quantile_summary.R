args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Read data (assuming the required column is the 5th)
data <- read.table(input_file, header = FALSE)
column_5 <- data$V5

# Compute the required quantiles
quantiles <- quantile(column_5, probs = c(0.95, 0.96, 0.97, 0.98, 0.99))

# Count the number of values in each quantile range
counts <- sapply(quantiles, function(q) sum(column_5 >= q))

# Extract population code from filename
population_code <- unlist(strsplit(basename(input_file), "_"))[1]

# Prepare the output
output <- data.frame(
  Population = population_code,
  TotalValues = length(column_5),
  Quantile95 = quantiles[1],
  Count95 = counts[1],
  Quantile96 = quantiles[2],
  Count96 = counts[2],
  Quantile97 = quantiles[3],
  Count97 = counts[3],
  Quantile98 = quantiles[4],
  Count98 = counts[4],
  Quantile99 = quantiles[5],
  Count99 = counts[5]
)

# Write the summary to the output file
write.table(output, file = output_file, row.names = FALSE, sep = "\t", quote = FALSE)
