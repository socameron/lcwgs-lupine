#!/usr/bin/env Rscript
# Usage: Rscript convert_vcf_to_genepop.R <input.vcf(.gz)> <output.gen>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript convert_vcf_to_genepop.R <input.vcf(.gz)> <output.gen>")
vcf_file    <- args[1]
genepop_out <- args[2]

# Threads (from SLURM if available)
thr_env <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")
threads <- suppressWarnings(as.integer(thr_env)); if (is.na(threads) || threads < 1) threads <- 1

message("Input VCF: ", vcf_file)
message("Output GENEPOP: ", genepop_out)
message("Threads: ", threads)

# Load radiator
suppressPackageStartupMessages({
  ok <- requireNamespace("radiator", quietly = TRUE)
})
if (!ok) stop("Package 'radiator' not installed. In R: install.packages('radiator')")

# Optional: strata file (INDIVIDUALS,STRATA) via env var
strata_fp <- Sys.getenv("RAD_STRATA", unset = "")
if (nzchar(strata_fp) && !file.exists(strata_fp)) {
  stop("RAD_STRATA set but file not found: ", strata_fp)
}

# Ensure output dir + unique work dir
outdir <- dirname(genepop_out)
if (!dir.exists(outdir) && !dir.create(outdir, recursive = TRUE)) {
  stop("Output directory is not writable: ", outdir)
}
run_dir <- tempfile("radiator_", tmpdir = outdir)
if (!dir.create(run_dir, recursive = TRUE, showWarnings = FALSE) && !dir.exists(run_dir)) {
  stop("Failed to create run directory: ", run_dir)
}
message("radiator working directory: ", run_dir)

# Build args for radiator (GT-only; skips AD/DP/PL parsing)
args_list <- list(
  data                = vcf_file,
  output              = "genepop",
  filename            = "genepop_export",
  path.folder         = run_dir,
  parallel.core       = threads,
  verbose             = TRUE,
  vcf.metadata        = FALSE,  # <-- keep only GT
  filter.common.markers = FALSE # faster for single-POP inputs
  # filter.monomorphic = TRUE   # keep if you want only polymorphic loci
)
if (nzchar(strata_fp)) args_list$strata <- strata_fp

message("Running radiator::genomic_converter() ...")
res <- try(do.call(radiator::genomic_converter, args_list), silent = TRUE)
if (inherits(res, "try-error")) {
  stop("radiator::genomic_converter failed:\n", as.character(res))
}

# Locate the produced .gen and move to final path
gen_files <- list.files(run_dir, pattern = "\\.gen$", full.names = TRUE, recursive = TRUE)
if (!length(gen_files)) stop("No .gen file produced under: ", run_dir)
gen_src <- gen_files[order(file.info(gen_files)$mtime, decreasing = TRUE)][1]

if (!file.copy(gen_src, genepop_out, overwrite = TRUE)) {
  stop("Failed to write final output to: ", genepop_out)
}
message("Wrote Genepop: ", genepop_out)
message("Done.")
