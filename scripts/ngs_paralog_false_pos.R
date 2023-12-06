## R script to create a graph from "{population}_paralog_WGS.lr" data ##

options(scipen = 999) # Prevent any scientific notation printing

# Process LR files
args <- commandArgs(trailingOnly = TRUE)
print(args)
input_file <- args[1]
output_file_bp0 <- args[2]
output_file_bp1 <- args[3]


# Print file paths for debugging
cat("Input file:", input_file, "\n")
cat("Output file bp0:", output_file_bp0, "\n")
cat("Output file bp1:", output_file_bp1, "\n")

lr <- read.table(input_file, header=FALSE, sep="\t")
lr$pval <- 0.5 * pchisq(lr$V5, df=1, lower.tail=FALSE)
lr$pval.adj <- p.adjust(lr$pval, method="fdr")

# Setting proportion of kept SNPs at 50%, controlling for false positives
# In other words, paralog.sites will be paralogous or deviant SNPs 
# Here we keep all 5 columns for bp.1 that we can input the .lr file for dupHMM (need all 5 columns - 1-indexed).
# Here we only keep 2 columns for bp.0 so we can create a BED file (0-indexed)
paralog.sites.bp0 <- lr[-which(lr$pval.adj < 0.1), 1:2]
paralog.sites.bp1 <- lr[-which(lr$pval.adj < 0.1), 1:5]


# Converting paralog.sites to BED format for filtering in SAMtools (0-based start)
# Note that the output of ngsParalog is 1-indexed
paralog.sites.bp0$V3 <- paralog.sites.bp0$V2
paralog.sites.bp0$V2 <- paralog.sites.bp0$V2 - 1

write.table(paralog.sites.bp0, file=output_file_bp0, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
write.table(paralog.sites.bp1, file=output_file_bp1, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)