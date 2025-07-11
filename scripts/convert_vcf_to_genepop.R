#!/usr/bin/env Rscript
# args: input.vcf output.gen
args <- commandArgs(trailingOnly=TRUE)
vcf_file    <- args[1]
genepop_out <- args[2]

library(vcfR)      # CRAN
library(adegenet)  # CRAN

# 1) Read VCF
vcf     <- read.vcfR(vcf_file)
# 2) Convert to genind
genind  <- vcfR2genind(vcf)
# 3) Write GENEPOP
#    by default, all samples are one “population”—you
#    can split them afterwards if you have multiple pops
genepop(genind, file = genepop_out)
