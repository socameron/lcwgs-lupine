#!/usr/bin/env Rscript
# estimate_RZooRoH.R
#
# Usage:
#   Rscript estimate_RZooRoH.R \
#     <input_zoo_file> \
#     <sample_file> \
#     <allele_freq_file> \
#     <output_base> \
#     <population> \
#     <n_threads>
#
# Example:
# Rscript estimate_RZooRoH.R \
#   results/.../APB_Zoo_format.txt \
#   data/lists/hap2/APB_sample_names.txt \
#   results/.../APB_allelefreq.txt \
#   results/.../APB/APB_summary \
#   APB \
#   4

suppressPackageStartupMessages({
  library(RZooRoH)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 6) {
  stop("Usage: Rscript estimate_RZooRoH.R <input_zoo> <sample_file> <allele_freq> <output_base> <population> <n_threads>")
}
input_zoo   <- args[1]
sample_file <- args[2]
allele_freq <- args[3]
output_base <- args[4]
population  <- args[5]
n_threads   <- as.integer(args[6])

output_dir <- dirname(output_base)
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 1) Load the data
zoo_data <- zoodata(file1      = input_zoo,
                    zformat    = "gl",
                    samplefile = sample_file,
                    allelefreq = allele_freq)

# 2) Define exactly the MixKR rate‐sets you want to try:
custom_models <- list(
  MixKR_3a = c(2,   16,    128),
  MixKR_3b = c(64,  256,   1024),
  MixKR_4a = c(8,   32,    128,  512),
  MixKR_4b = c(128, 256,   512,  1024),
  MixKR_5a = c(8,   16,    32,   64,   128),
  MixKR_5b = c(128, 256,   512,   1024, 2048),
  MixKR_5c = c(4,   16,    512,  2048, 8192),
  MixKR_6a = c(8,   32,    512,  2048, 8192,  32768),
  MixKR_6b = c(64,  128,   256,  512,  1024,  2048)
)

candidate_models <- setNames(
  lapply(custom_models, function(rates){
    zoomodel(predefined = TRUE,
             K          = length(rates),
             krates     = rates,
             layers     = FALSE)
  }),
  names(custom_models)
)

# 3) Fit each, record total BIC, pick best
bic_table <- data.frame(
  Population = character(0),
  Model      = character(0),
  K          = integer(0),
  Rates      = character(0),
  TotalBIC   = numeric(0),
  stringsAsFactors = FALSE
)
best_bic    <- Inf
best_name   <- NULL
best_result <- NULL

for(mn in names(candidate_models)){
  cat("→ Running model:", mn, "…\n")
  zmod <- candidate_models[[mn]]
  res  <- tryCatch(
    zoorun(zmod, zoo_data, nT = n_threads),
    error = function(e){
      warning("  Model ", mn, " failed: ", e$message)
      NULL
    }
  )
  if(is.null(res)) next

  tbic <- sum(res@modbic)
  bic_table <- rbind(bic_table, data.frame(
    Population = population,
    Model      = mn,
    K          = zmod@K,
    Rates      = paste(zmod@krates, collapse=","),
    TotalBIC   = tbic,
    stringsAsFactors = FALSE
  ))

  if(tbic < best_bic){
    best_bic    <- tbic
    best_name   <- mn
    best_result <- res
  }
}

# 4) Save model comparison
cmp_file <- file.path(output_dir, paste0(population, "_model_comparison.csv"))
write.csv(bic_table, cmp_file, row.names = FALSE)
cat("✔ Model comparison →", cmp_file, "\n")

# 5) Export best‐model results
if(is.null(best_result)) stop("❌ No model converged successfully!")

# 5a) realized autozygosity per class
realized_df <- as.data.frame(realized(best_result))
real_file   <- file.path(output_dir, paste0(population, "_realized_", best_name, ".csv"))
write.csv(realized_df, real_file, row.names = FALSE)
cat("✔ Realized autozygosity →", real_file, "\n")

# 5b) genome‐wide inbreeding F_{G–T} at T = max rate
maxR    <- max(best_result@krates)
fg      <- cumhbd(best_result, maxR)
fg_file <- file.path(output_dir, paste0(population, "_inbreeding_", best_name, ".csv"))
write.csv(data.frame(F = fg), fg_file, row.names = FALSE)
cat("✔ Inbreeding coefficients →", fg_file, "\n")

# 5c) HBD segments
segs     <- best_result@hbdseg
segs_file <- file.path(output_dir, paste0(population, "_HBDsegments_", best_name, ".csv"))
write.csv(segs, segs_file, row.names = FALSE)
cat("✔ HBD segments →", segs_file, "\n\n")

# 6) Plots — dump into a single multi‐page PDF
plots_pdf <- file.path(output_dir, paste0(population, "_diagnostics.pdf"))
pdf(plots_pdf, width=8, height=6)

# 6.1 Population-level barplot of F_G by class
zooplot_prophbd(list(POP = best_result),
                cols  = 'tomato',
                style = 'barplot')

# 6.2 Population-level cumulative inbreeding curve
zooplot_prophbd(list(POP = best_result),
                style      = 'lines',
                cumulative = TRUE)

# 6.3 Individual-level cumulative inbreeding
zooplot_individuals(list(POP = best_result),
                    cumulative = TRUE)

# 6.4 Individual‐level partitioning (stacked bar)
zooplot_partitioning(list(POP = best_result),
                    ylim    = c(0, 0.5),
                    nonhbd  = FALSE,
                    plotids = FALSE)

# 6.5 A few HBD segments on “Scaffold_1” from 1–2 Mb
# (just as an example, adjust to your real scaffold names & region)
zooplot_hbdseg(list(POP = best_result),
               randomids = TRUE,
               nrandom   = c(5),
               chr       = 1,
               coord     = c(1e6, 2e6))

dev.off()
cat("✔ Diagnostic plots →", plots_pdf, "\n")
cat("Done for", population, "— best model was", best_name, "\n")
