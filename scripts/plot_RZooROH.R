#!/usr/bin/env Rscript
#
# plot_RZooROH.R
#
# USAGE:
#   Rscript plot_RZooROH.R \
#     <real1.csv> [<real2…>] \
#     <meta.csv> \
#     <recent_thr> \
#     <out.pdf> \
#     <png_prefix>
#

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 6) {
  stop("Usage: plot_RZooROH.R <real1>… <meta> <recent_thr> <out.pdf> <png_prefix>")
}

# parse
png_prefix   <- tail(args,1)
out_pdf      <- tail(args,2)[1]
recent_thr   <- as.numeric(tail(args,3)[1])
meta_fp      <- tail(args,4)[1]
realized_fps <- head(args, -4)

# your colours & theme
cols <- c(Core="#9C67BD", Edge="#43AEBD")
population_colors <- c(
  "C.OH.1" = "#DA70D6",    "C.MI.2" = "#9932CC",
  "C.ON.3" = "#551A8B",    "C.IN.2" = "#B19CD9",
  "C.IL.1" = "#9370DB",    "C.MN.4" = "#BA55D3",
  "C.WI.1" = "#663399",    "C.WI.2" = "#9400D3",
  "E.ON.5" = "#00FFFF",    "E.ON.8" = "#7FFFD4",
  "E.ON.9" = "#00868B",    "E.MI.7" = "#20B2AA",
  "E.MN.2" = "#48D1CC",    "E.MN.3" = "#5F9EA0",
  "E.NH.1" = "#008B8B",    "E.NY.1" = "#00CED1",
  "E.WI.4" = "#40E0D0",    "E.WI.9" = "#009999"
)
custom_theme <- theme(
  panel.border      = element_blank(),
  panel.grid.major  = element_blank(),
  panel.grid.minor  = element_blank(),
  axis.line         = element_line(colour="black", linewidth=1),
  axis.ticks        = element_line(linewidth=1),
  axis.ticks.length = unit(0.25, "cm"),
  plot.title        = element_text(hjust=0.5),
  axis.text         = element_text(size=14, angle=45, hjust=1),
  axis.title        = element_text(size=14),
  legend.position   = "none"
)
theme_set(custom_theme)

# 1) metadata
meta <- read.csv(meta_fp, stringsAsFactors=FALSE) %>%
  select(Site_code, Site_code2, Range_position) %>%
  mutate(
    Label          = sub("^[CE]\\.","", Site_code2),
    Range_position = factor(Range_position, c("Core","Edge"))
  )

# 2) read realized & compute F_recent / F_ancient
df_list <- lapply(realized_fps, function(fp){
  pop <- sub("_realized_.*\\.csv$","", basename(fp))
  dat <- read.csv(fp, check.names=FALSE)
  hbd <- grep("^R_[0-9]+$", names(dat), value=TRUE)
  rates <- as.numeric(sub("^R_","", hbd))
  rec   <- hbd[rates <= recent_thr]
  anc   <- hbd[rates  > recent_thr]
  data.frame(
    Site_code = pop,
    F_recent  = rowSums(dat[ , rec, drop=FALSE]),
    F_ancient = rowSums(dat[ , anc, drop=FALSE]),
    stringsAsFactors=FALSE
  )
})
plot_df <- bind_rows(df_list) %>%
  inner_join(meta, by="Site_code")

# 3) build plots
p1 <- ggplot(plot_df, aes(Range_position, F_recent, fill=Range_position)) +
  geom_boxplot() +
  scale_fill_manual(values=cols) +
  labs(y="Recent autozygosity (F_recent)", x=NULL) +
  custom_theme

p2 <- ggplot(plot_df, aes(Range_position, F_ancient, fill=Range_position)) +
  geom_boxplot() +
  scale_fill_manual(values=cols) +
  labs(y="Ancient autozygosity (F_ancient)", x=NULL) +
  custom_theme

pop_sum <- plot_df %>%
  group_by(Site_code, Label, Range_position) %>%
  summarize(
    mFrec = mean(F_recent),
    mFanc = mean(F_ancient),
    .groups="drop"
  )

p3 <- ggplot(pop_sum, aes(mFanc, mFrec, color=Range_position)) +
  geom_point(size=3) +
  geom_text_repel(
    aes(label = Label),
    size         = 3,
    box.padding  = 0.4,
    point.padding= 0.2,
    segment.size = 0.2
  ) +
  scale_color_manual(values=cols) +
  labs(x = expression(paste("Mean ancient ", F[ROH])), y = expression(paste("Mean recent ", F[ROH]))) +
  custom_theme
  

# 4) write PDF
pdf(out_pdf, width=7, height=6)
print(p1); print(p2); print(p3)
dev.off()

# 5) write three independent PNGs
ggsave(paste0(png_prefix, "_Frecent.png"), p1, width=6, height=5, dpi=150)
ggsave(paste0(png_prefix, "_Fancient.png"), p2, width=6, height=5, dpi=150)
ggsave(paste0(png_prefix, "_scatter.png"),  p3, width=6, height=5, dpi=150)

message("✅ Written:\n  PDF: ", out_pdf,
        "\n  PNGs:\n    ", paste0(png_prefix, "_Frecent.png"),
        "\n    ", paste0(png_prefix, "_Fancient.png"),
        "\n    ", paste0(png_prefix, "_scatter.png"))

