## R script to create a graph from "{population}_paralog_WGS.lr" data ##

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

# Setting proportion of kept SNPs at 5 or 10%, controlling for false positives
# In other words, paralog.sites will be paralogous or deviant SNPs
paralog.sites <- lr[-which(lr$pval.adj < 0.10), 1:2]


# Converting paralog.sites to BED format for filtering in SAMtools (0-based start)
paralog.sites.bp0 <- paralog.sites
paralog.sites.bp0$V3 <- paralog.sites$V2
paralog.sites.bp0$V2 <- paralog.sites.bp0$V2 - 1

write.table(paralog.sites.bp0, file=output_file_bp0, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
write.table(paralog.sites, file=output_file_bp1, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)