#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(ggrepel)
})

# --------- CLI ---------
parser <- ArgumentParser()
parser$add_argument("--nc_csv", type="character", required=TRUE)
parser$add_argument("--nsites_q", type="double", default=0.15)
parser$add_argument("--nc_year", type="character", required=TRUE) # "2023"/"2024"/"2025"
parser$add_argument("--scatter", type="character", required=TRUE)
parser$add_argument("--summary", type="character", required=TRUE)
args <- parser$parse_args()

nc_csv   <- args$nc_csv
ns_q     <- args$nsites_q
nc_year  <- args$nc_year
out_plot <- args$scatter
out_sum  <- args$summary

# --------- IDs & mapping (old -> new 'Site_code2') ---------
old_populations <- c("HPW","IDNP-MW","LCTGP","MFNP","PPP","RLPLV","SWCP",
                     "APB","BSL","BSNA","CPB","FMB","GRAY","NBWA","NGP","PBBT","RSB","UWA")

pop_mapping <- c(
  "NGP"     = "C.IL.1",
  "IDNP-MW" = "C.IN.2",
  "MFNP"    = "C.MI.2",
  "RSB"     = "C.MN.4",
  "LCTGP"   = "C.OH.1",
  "SWCP"    = "C.ON.3",
  "UWA"     = "C.WI.1",
  "FMB"     = "C.WI.2",
  "GRAY"    = "E.MI.7",
  "PBBT"    = "E.MN.2",
  "BSL"     = "E.MN.3",
  "CPB"     = "E.NH.1",
  "APB"     = "E.NY.1",
  "HPW"     = "E.ON.5",
  "RLPLV"   = "E.ON.8",
  "PPP"     = "E.ON.9",
  "BSNA"    = "E.WI.4",
  "NBWA"    = "E.WI.9"
)

population_colors <- c(
  "C.OH.1"="#DA70D6","C.MI.2"="#9932CC","C.ON.3"="#551A8B","C.IN.2"="#B19CD9",
  "C.IL.1"="#9370DB","C.MN.4"="#BA55D3","C.WI.1"="#663399","C.WI.2"="#9400D3",
  "E.ON.5"="#00FFFF","E.ON.8"="#7FFFD4","E.ON.9"="#00868B","E.MI.7"="#20B2AA",
  "E.MN.2"="#48D1CC","E.MN.3"="#5F9EA0","E.NH.1"="#008B8B","E.NY.1"="#00CED1",
  "E.WI.4"="#40E0D0","E.WI.9"="#009999"
)

# --------- helpers ---------
read_theta_pestPG <- function(file_path) {
  df <- suppressWarnings(read.table(file_path, header = TRUE, sep = "\t",
                                    comment.char = "#", skip = 1,
                                    stringsAsFactors = FALSE))
  colnames(df) <- c("Region","Chr","WinCenter","tW","tP","tF","tH","tL",
                    "TajimaD","FuLiF","FuLiD","FayH","Zeng","nSites")
  df %>%
    transmute(
      window_position = WinCenter,
      tajimasD = TajimaD,
      tP = tP,
      nSites = nSites
    )
}

pi_per_window <- function(tp, nsites) {
  out <- ifelse(nsites > 0, tp / nsites, NA_real_)
  as.numeric(out)
}

# --------- ingest Nc table ---------
nc_df <- read_csv(nc_csv, show_col_types = FALSE)
# Expect: Site_code, Site_code2, Range_position, 2023_Nc, 2024_Nc, 2025_Nc
nc_col <- paste0(nc_year, "_Nc")
if (!nc_col %in% names(nc_df)) stop("Nc column ", nc_col, " not found in ", nc_csv)

# --------- read all hap2 files ---------
combined <- list()

for (pop_old in old_populations) {
  new_pop <- unname(pop_mapping[pop_old])
  file_path <- file.path("results/theta/hap2", sprintf("%s_out.thetasWindow.gz.pestPG", pop_old))

  if (!file.exists(file_path)) {
    warning("Missing: ", file_path)
    next
  }

  th <- read_theta_pestPG(file_path)

  # per-pop window filter by nSites
  qcut <- suppressWarnings(quantile(th$nSites, ns_q, na.rm = TRUE))
  thf <- th %>% filter(is.finite(nSites), nSites > qcut)

  thf <- thf %>%
    mutate(
      population = new_pop,
      pi_window  = pi_per_window(tP, nSites)
    )

  combined[[length(combined) + 1]] <- thf
}

if (length(combined) == 0) stop("No theta data loaded from hap2.")
theta_all <- bind_rows(combined)

# --------- summarise π by population ---------
pop_summary <- theta_all %>%
  group_by(population) %>%
  summarise(
    n_windows = sum(!is.na(pi_window)),
    pi_mean   = mean(pi_window, na.rm = TRUE),
    pi_median = median(pi_window, na.rm = TRUE),
    pi_sd     = sd(pi_window, na.rm = TRUE),
    .groups = "drop"
  )

# --------- join to Nc (by Site_code2) ---------
joined <- pop_summary %>%
  left_join(
    nc_df %>% select(Site_code2, Range_position, all_of(nc_col)),
    by = c("population" = "Site_code2")
  ) %>%
  rename(Nc = all_of(nc_col)) %>%          # <— renames "2025_Nc" => "Nc"
  mutate(group = ifelse(startsWith(population, "C."), "Core", "Edge"))

# write summary csv
readr::write_csv(joined, out_sum)

# --------- plot π vs Nc ---------
plot_df <- joined %>%
  mutate(Nc = suppressWarnings(as.numeric(Nc)),
         pi_mean = as.numeric(pi_mean)) %>%
  filter(is.finite(Nc), is.finite(pi_mean))

# data to FIT (exclude the huge pop), but still SHOW all points
fit_df <- dplyr::filter(plot_df, population != "C.IN.2")

# fit and equation using fit_df
fit <- lm(pi_mean ~ Nc, data = fit_df)
eq_label <- sprintf("y = %.6g + %.6g · x,   R² = %.3f",
                    coef(fit)[1], coef(fit)[2], summary(fit)$r.squared)

# x-limits based on non-outlier pops
xlims <- range(fit_df$Nc, na.rm = TRUE)

p <- ggplot(plot_df, aes(x = Nc, y = pi_mean)) +
  # all points (including C.IN.2) in pop colors
  geom_point(aes(color = population), size = 3) +
  # regression line/SE from fit_df only
  geom_smooth(data = fit_df, aes(x = Nc, y = pi_mean),
              method = "lm", se = TRUE, color = "black", linewidth = 0.9) +
  # labels same color as points, no leader lines
  ggrepel::geom_text_repel(
    aes(label = population, color = population),
    size = 3, segment.color = NA, max.overlaps = Inf, box.padding = 0.25, seed = 123
  ) +
  scale_color_manual(values = population_colors, guide = "none") +
  coord_cartesian(xlim = xlims) +   # zoom to non-outlier range
  labs(
    title = paste0("Nucleotide diversity (", "\u03C0", ") vs Nc (", nc_year, ")"),
    x = paste0("Nc (", nc_year, ")"),
    y = expression("Mean nucleotide diversity " * pi)
  ) +
  theme_classic(base_size = 13) +
  theme(axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 0.8)) +
  annotate("text", x = Inf, y = -Inf, label = eq_label,
           hjust = 1.02, vjust = -0.6)
           
ggsave(filename = out_plot, plot = p, width = 7.5, height = 5.5, dpi = 300)
