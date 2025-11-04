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
  make_option("--n_runs", type="integer",   default=10L,
              help="Number of multi-start runs for MLE [default %default]."),
  make_option("--models", type="character", default="gammaexp",
              help="Comma-separated: gammaexp,discrete,gammadiscrete,displacedgamma or 'all'. [default %default]"),
  make_option("--bootstrap", type="logical", default=FALSE,
              help="If TRUE, enable bootstrapping to get CIs on discretized masses [default %default]."),
  make_option("--ci_level", type="double", default=0.05,
              help="Alpha level for CIs (e.g., 0.05 -> 95%% CI) [default %default]."),
  make_option("--bootstrap_type", type="character", default="percentile",
              help="percentile or bca [default %default]"),
  make_option("--fix_beneficial_for_folded", type="logical", default=TRUE,
              help="Fix beneficial part (p_b=0, S_b=1, eps=0) for folded SFS [default %default].")
)

opt <- parse_args(OptionParser(option_list = opt_list))
if (is.null(opt$neut) || is.null(opt$sel) || is.null(opt$outdir) || is.null(opt$pop)) {
  stop("Missing required args: --neut, --sel, --outdir, --pop")
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
png_fp      <- file.path(opt$outdir, sprintf("%s_fastdfe_models_plot.png", opt$pop))
json_fp     <- file.path(opt$outdir, sprintf("%s_fastdfe_inference.json", opt$pop))
csv_fp      <- file.path(opt$outdir, sprintf("%s_fastdfe_discretized_models.csv", opt$pop))
summary_txt <- file.path(opt$outdir, sprintf("%s_fastdfe_summary.txt", opt$pop))

# ---------- Load fastDFE ----------
fd <- load_fastdfe()
BaseInference <- fd$BaseInference
Spectrum      <- fd$Spectrum
InferenceAgg  <- fd$Inference   # aggregator for plotting / summarizing across inferences (R wrapper)

# ---------- helpers ----------
read_sfs <- function(path) {
  if (!file.exists(path)) stop(sprintf("File not found: %s", path))
  x <- scan(path, what = numeric(), quiet = TRUE)
  x[is.finite(x)]
}

# map CLI tokens -> constructor functions
make_parametrizations <- function(fd, models_arg) {
  key <- tolower(trimws(models_arg))
  if (identical(key, "all")) key <- "gammaexp,discrete,gammadiscrete,displacedgamma"
  toks <- unique(strsplit(key, ",", fixed = TRUE)[[1]])
  toks <- trimws(tolower(toks))
  ctor <- list(
    gammaexp       = fd$GammaExpParametrization,
    discrete       = fd$DiscreteParametrization,
    gammadiscrete  = fd$GammaDiscreteParametrization,
    displacedgamma = fd$DisplacedGammaParametrization
  )
  unknown <- setdiff(toks, names(ctor))
  if (length(unknown)) stop("Unknown model(s): ", paste(unknown, collapse = ", "))
  setNames(lapply(toks, function(t) ctor[[t]]()), toks)
}

# discretization bins (same idea as your script; add +Inf so we can call the aggregator if needed)
intervals_num <- c(-Inf, -100, -10, -1, 0, Inf)

# For folded SFS, it is common to fix beneficial part; keep your defaults but make it switchable
fixed <- if (isTRUE(opt$fix_beneficial_for_folded)) {
  list(all = list(eps = 0, S_b = 1, p_b = 0))
} else NULL

# ---------- Data ----------
sfs_neut <- Spectrum(c(read_sfs(opt$neut)))
sfs_sel  <- Spectrum(c(read_sfs(opt$sel)))

# ---------- Build and run one or many models ----------
params_list <- make_parametrizations(fd, opt$models)
labels <- names(params_list)

inferences <- vector("list", length(labels))
names(inferences) <- labels

for (lab in labels) {
  inf <- BaseInference(
    sfs_neut     = sfs_neut,
    sfs_sel      = sfs_sel,
    n_runs       = as.integer(opt$n_runs),
    model        = params_list[[lab]],
    do_bootstrap = isTRUE(opt$bootstrap),
    fixed_params = fixed
  )
  inf$run()
  inferences[[lab]] <- inf
}

# ---------- Plot (multi-model) ----------
# Use the R wrapper's aggregator; same as docs example
# Note: plot_discretized() handles multiple inferences and labels
p <- InferenceAgg$plot_discretized(
  inferences = unname(inferences),
  labels     = labels
)

if (inherits(p, "gg") || inherits(p, "ggplot")) {
  ggsave(filename = png_fp, plot = p, width = 7.2, height = 4.6, dpi = 220)
} else {
  # matplotlib fallback
  plt <- reticulate::import("matplotlib.pyplot", convert = TRUE)
  fig <- tryCatch({
    if (reticulate::py_has_attr(p, "savefig")) p else plt$gcf()
  }, error = function(e) NULL)
  if (!is.null(fig) && reticulate::py_has_attr(fig, "savefig")) {
    fig$savefig(png_fp, dpi = 220L, bbox_inches = "tight")
  }
}

# ---------- Discretized masses table (optionally with CIs if bootstrapped) ----------
# If bootstrapping is enabled, use aggregator to fetch values + errors; else call per inference
get_discretized_safe <- function() {
  if (!isTRUE(opt$bootstrap)) {
    # no CIs: collect per inference
    out <- do.call(rbind, lapply(names(inferences), function(lab) {
      inf <- inferences[[lab]]
      masses <- tryCatch(as.numeric(inf$get_discretized(intervals = intervals_num)), error = function(e) NULL)
      if (is.null(masses)) masses <- rep(NA_real_, length(intervals_num) - 1L)
      data.frame(
        pop       = opt$pop,
        model     = lab,
        bin_left  = head(intervals_num, -1),
        bin_right = tail(intervals_num, -1),
        mass      = masses,
        ci_low    = NA_real_,
        ci_high   = NA_real_
      )
    }))
    return(out)
  } else {
    # use aggregator to compute values and CI across inferences
    # R wrapper mirrors Python static method signature
    res <- InferenceAgg$get_discretized(
      inferences = unname(inferences),
      labels     = labels,
      intervals  = intervals_num,
      confidence_intervals = TRUE,
      ci_level   = as.numeric(opt$ci_level),
      bootstrap_type = opt$bootstrap_type
    )
    # 'res' typically returns a dict-like object with $values and $errors (or similar)
    # Harmonize into a data.frame defensively
    vals <- res$values
    errs <- res$errors
    do.call(rbind, lapply(labels, function(lab) {
      v <- as.numeric(vals[[lab]])
      e <- if (!is.null(errs)) as.numeric(errs[[lab]]) else rep(NA_real_, length(v))
      # Convert symmetric errors to CI bounds around v (approx)
      df <- data.frame(
        pop       = opt$pop,
        model     = lab,
        bin_left  = head(intervals_num, -1),
        bin_right = tail(intervals_num, -1),
        mass      = v,
        ci_low    = ifelse(is.finite(e), v - e, NA_real_),
        ci_high   = ifelse(is.finite(e), v + e, NA_real_)
      )
      df
    }))
  }
}

disc_df <- get_discretized_safe()
write.csv(disc_df, csv_fp, row.names = FALSE)

# ---------- JSON metadata ----------
meta <- list(
  pop       = opt$pop,
  n_runs    = opt$n_runs,
  models    = labels,
  intervals = intervals_num,
  bootstrap = opt$bootstrap,
  ci_level  = opt$ci_level,
  bootstrap_type = opt$bootstrap_type,
  fixed_params  = if (is.null(fixed)) NULL else fixed
)
cat(jsonlite::toJSON(meta, auto_unbox = TRUE, null = "null", pretty = TRUE),
    file = json_fp)

# ---------- Human-readable summary ----------
sink(summary_txt)
cat("fastDFE multi-model summary\n")
cat("Population:", opt$pop, "\n")
cat("n_runs:", opt$n_runs, "\n")
cat("Models:", paste(labels, collapse = ", "), "\n")
cat("Intervals:", paste(intervals_num, collapse = ", "), "\n")
cat("Bootstrap:", opt$bootstrap, " type:", opt$bootstrap_type, " alpha:", opt$ci_level, "\n\n")
for (lab in labels) {
  cat("---- Model:", lab, "----\n")
  print(inferences[[lab]])
  cat("\n")
}
sink()

message("Done: ", opt$pop)
