args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Read the .lr file (adjust read.table as necessary based on your file format)
lr_data <- read.table(input_file, header = FALSE, colClasses = c("character", "numeric"))

# Assuming the first column is chromosome (chr) and the second column is position (pos)
# BED format: chr, start (0-based), end (exclusive)
lr_data$V2 <- lr_data$V2 - 1 # Convert to 0-based indexing for start position

# Use formatC to avoid scientific notation
bed_data <- data.frame(chr = lr_data$V1, 
                       start = formatC(lr_data$V2, format = "f", digits = 0), 
                       end = formatC(lr_data$V2 + 1, format = "f", digits = 0))

# Write the BED file
write.table(bed_data, 
            file = output_file, 
            row.names = FALSE, 
            col.names = FALSE, 
            sep = "\t", 
            quote = FALSE,
            eol = "\n")
