#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
out_plot  <- tail(args, 1)
in_files  <- head(args, -1)

library(dplyr); library(ggplot2); library(readr)

# *** PARAMETERS you choose ***
Ne              <- 5000             # your assumed effective population size
recomb_per_bp   <- 1e-6             # 1 cM/Mb → 1e-8 per bp, adjust if needed
sample_size     <- 25               # number of diploid individuals
# Hill & Weir finite sample correction factor: 1/(n) for diploids use 1/(2N)? 
# (approximate, or use the full formula if you have it)

# Read & combine samples
ld_df <- bind_rows(lapply(in_files, function(f){
  pop <- sub("_ld_pair_sample.txt$", "", basename(f))
  df  <- read_tsv(f, col_names=c("dist_bp","r2"), show_col_types=FALSE)
  df$population <- pop
  df
}))

# Convert distances to Mb
ld_df <- ld_df %>% mutate(dist_Mb = dist_bp / 1e6)

# Bin by distance (e.g. 100 windows across range)
ld_binned <- ld_df %>%
  group_by(population, bin = cut(dist_Mb, breaks=seq(0, max(dist_Mb), length=100))) %>%
  summarize(
    mid = mean(as.numeric(sub("\\((.+),(.+)\\]", "\\1", bin)) + 
               diff(as.numeric(strsplit(as.character(bin),",")[[1]]))/2),
    mean_r2 = mean(r2, na.rm=TRUE)
  )

# Theoretical curve: (10 + C)/((2 + C)*(11 + C))
Cfunc <- function(d) 4 * Ne * d * recomb_per_bp
expect_r2 <- function(d) {
  C <- Cfunc(d)
  (10 + C) / ((2 + C)*(11 + C))
}

# Plot
p <- ggplot() +
  geom_point(data=ld_binned, aes(x=mid, y=mean_r2), alpha=0.6) +
  facet_wrap(~population, scales="free_y", ncol=4) +
  stat_function(fun=expect_r2, aes(colour="Theory"), size=1) +
  labs(x="Distance (Mb)", y=expression(mean~r^2),
       title="Observed vs. Hill & Weir expectation") +
  theme_minimal() +
  theme(legend.position="bottom")

ggsave(out_plot, p, width=12, height=8)
