#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(fastdfe)   # R wrapper over Python
  library(ggplot2)   # for ggsave
  library(reticulate)
})

# ---------- CLI ----------
opt_list <- list(
  make_option("--neut",   type="character", help="Path to folded neutral SFS (4-fold)."),
  make_option("--sel",    type="character", help="Path to folded selected SFS (CDS)."),
  make_option("--outdir", type="character", help="Output directory."),
  make_option("--pop",    type="character", help="Population code."),
  make_option("--n_runs", type="integer",   default=10L)
)
opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$neut) || is.null(opt$sel) || is.null(opt$outdir) || is.null(opt$pop)) {
  stop("Missing required args: --neut, --sel, --outdir, --pop")
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
png_fp      <- file.path(opt$outdir, sprintf("%s_fastdfe_plot.png", opt$pop))
json_fp     <- file.path(opt$outdir, sprintf("%s_fastdfe_inference.json", opt$pop))
csv_fp      <- file.path(opt$outdir, sprintf("%s_fastdfe_discretized.csv", opt$pop))
summary_txt <- file.path(opt$outdir, sprintf("%s_fastdfe_summary.txt", opt$pop))

# ---------- Load fastDFE classes via the R wrapper ----------
fd <- load_fastdfe()
BaseInference <- fd$BaseInference
Spectrum      <- fd$Spectrum

# ---------- helpers ----------
read_sfs <- function(path) {
  if (!file.exists(path)) stop(sprintf("File not found: %s", path))
  x <- scan(path, what = numeric(), quiet = TRUE)
  x[is.finite(x)]
}

# ---------- Data ----------
sfs_neut <- Spectrum(c(read_sfs(opt$neut)))
sfs_sel  <- Spectrum(c(read_sfs(opt$sel)))

# Built-in default bins for discretized plot/table
intervals_num <- c(-Inf, -100, -10, -1, 0)

# For folded SFS, fix beneficial component + epsilon (common practice)
fixed <- list(all = list(eps = 0, S_b = 1, p_b = 0))

# ---------- Inference ----------
inf <- BaseInference(
  sfs_neut     = sfs_neut,
  sfs_sel      = sfs_sel,
  n_runs       = as.integer(opt$n_runs),
  fixed_params = fixed
)
BaseInference$run(inf)

# ---------- Plot (ggplot in R wrapper; fallback to mpl if needed) ----------
p <- BaseInference$plot_discretized(inf, intervals = intervals_num)

if (inherits(p, "gg") || inherits(p, "ggplot")) {
  ggsave(filename = png_fp, plot = p, width = 6.5, height = 4.5, dpi = 200)
} else {
  plt <- reticulate::import("matplotlib.pyplot", convert = TRUE)
  fig <- tryCatch({
    if (reticulate::py_has_attr(p, "savefig")) p else plt$gcf()
  }, error = function(e) NULL)
  if (!is.null(fig) && reticulate::py_has_attr(fig, "savefig")) {
    fig$savefig(png_fp, dpi = 200L, bbox_inches = "tight")
  }
}

# ---------- Discretized masses table ----------
masses <- NULL
if (!is.null(inf$get_discretized)) {
  masses <- tryCatch(as.numeric(inf$get_discretized(intervals = intervals_num)), error = function(e) NULL)
}
if (is.null(masses)) masses <- rep(NA_real_, length(intervals_num) - 1L)

write.csv(
  data.frame(
    pop = opt$pop,
    bin_left  = head(intervals_num, -1),
    bin_right = tail(intervals_num, -1),
    mass      = masses
  ),
  csv_fp, row.names = FALSE
)


# ---------- JSON metadata (safe: no Python objects) ----------
meta <- list(
  pop       = opt$pop,
  n_runs    = opt$n_runs,
  intervals = intervals_num
)
cat(jsonlite::toJSON(meta, auto_unbox = TRUE, null = "null", pretty = TRUE),
    file = json_fp)


# ---------- Human-readable summary ----------
sink(summary_txt)
cat("fastDFE BaseInference summary\n")
cat("Population:", opt$pop, "\n")
cat("n_runs:", opt$n_runs, "\n")
cat("Intervals:", paste(intervals_num, collapse = ", "), "\n\n")
print(inf)
sink()

message("Done: ", opt$pop)
