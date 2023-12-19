# filter_gff.R
args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
output_with_col3 <- args[2]
output_without_col3 <- args[3]

options(scipen = 999)

gff_data <- read.table(input_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

cds_data <- gff_data[gff_data$V3 == "CDS", ]

# Extracting scaffold numbers and converting them to numeric
cds_data$scaffold_num <- as.numeric(gsub("Scaffold_([0-9]+).*", "\\1", cds_data$V1))

# Sorting by scaffold number and then by start position
cds_data <- cds_data[order(cds_data$scaffold_num, cds_data$V4), ]

# Preparing data for output
cds_with_col3 <- cds_data[, c(1, 3, 4, 5)]
cds_without_col3 <- cds_data[, c(1, 4, 5)]

write.table(cds_with_col3, file = output_with_col3, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, eol = "\n")
write.table(cds_without_col3, file = output_without_col3, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, eol = "\n")
