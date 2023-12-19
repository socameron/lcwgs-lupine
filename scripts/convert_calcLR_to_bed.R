args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]

# Read the .lr file without headers
lr_data <- read.table(input_file, header = FALSE)

# Assign column names
colnames(lr_data) <- c('chr', 'pos', 'neg_log_null', 'neg_log_alt', 'lr')

# Create BED format data
# BED format: chr, start (0-based), end (exclusive)
bed_data <- data.frame(
  chr = lr_data$chr,
  start = lr_data$pos - 1,  # Convert to 0-based indexing
  end = lr_data$pos         # End is start + 1, since BED format is half-open
)

# Write the BED file
write.table(bed_data, 
            file = output_file, 
            row.names = FALSE, 
            col.names = FALSE, 
            sep = "\t", 
            quote = FALSE,
            eol = "\n")

