#!/usr/bin/env Rscript
# args: input.vcf output.gen
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript vcf_to_genepop.R <input.vcf(.gz)> <output.genepop>")
}
vcf_file    <- args[1]
genepop_out <- args[2]


library(vcfR)      # CRAN
library(adegenet)  # CRAN

# 1) Read VCF
message("Reading VCF: ", vcf_file)
vcf     <- read.vcfR(vcf_file)
# 2) Convert to genind
message("Converting to genind...")
genind  <- vcfR2genind(vcf)
# 3) Write GENEPOP
#    by default, all samples are one “population”—you
#    can split them afterwards if you have multiple pops
message("Writing GENEPOP to: ", genepop_out)
genepop(genind, file = genepop_out)

message("Done.")