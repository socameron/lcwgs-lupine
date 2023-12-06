args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
BH_VAR <- as.numeric(args[3])

lr <- read.table(input_file, header=FALSE, sep="\t")
lr$pval <- 0.5 * pchisq(lr$V5, df=1, lower.tail=FALSE)
lr$pval.adj <- p.adjust(lr$pval, method="fdr")

# Adjusting the threshold based on BH_VAR
threshold <- BH_VAR / 100
paralog_sites <- lr[-which(lr$pval.adj < threshold), 1:2]

write.table(paralog_sites, file=output_file, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

