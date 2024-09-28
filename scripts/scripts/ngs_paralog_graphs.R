## R script to create a graph from "{population}_paralog_WGS.lr" data ##

# Process LR files
args <- commandArgs(trailingOnly = TRUE)
print(args)
input_file <- args[1]
plot_file <- args[2]

# Print file paths for debugging
cat("Input file:", input_file, "\n")
cat("Plot file:", plot_file, "\n")

lr <- read.table(input_file, header=FALSE, sep="\t")
lr$pval <- 0.5 * pchisq(lr$V5, df=1, lower.tail=FALSE)
lr$pval.adj <- p.adjust(lr$pval, method="fdr")

qc.sites <- lr[-which(lr$pval.adj < 0.05), 1:2]


# Load necessary libraries for Manhattan plot
library(ggplot2)

ggplot(lr, aes(x=V2, y=-log10(pval))) + geom_point() + theme_classic() + 
  labs(title="Manhattan Plot of Likelihood Ratios", x="Position", y="-log10(p-value)")
ggsave(plot_file, width = 10, height = 8)

