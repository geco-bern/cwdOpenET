# Plots the CONUS-wide CWDmax map together with CONUS-wide gap-filled months and
# a gap-filled months distribution (barplot)


# ---- Setup ------------------------------------------------------------------
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyterra)
library(scico)
library(scales)
library(here)


# Load data 
r_conus      <- terra::rast(here("data", "cwd_conus.tif"))
r_cwd        <- r_conus[["cwd_max"]]
r_gap_filled <- r_conus[["gap_filled_months"]]

# Downsample for plotting 
# 30 m -> 600 m (factor 20). Continuous variables -> mean aggregation.
r_cwd_agg        <- terra::aggregate(r_cwd,        fact = 20, fun = "mean")
r_gap_filled_agg <- terra::aggregate(r_gap_filled, fact = 20, fun = "mean")


# ---- Map: CWDmax ------------------------------------------------------------
p_map_cwd <- ggplot() +
  tidyterra::geom_spatraster(data = r_cwd_agg) +
  scale_fill_viridis_c(
    name      = expression(CWD[max]~"[mm]"),
    option    = "viridis",
    direction = -1,              # low = yellow, high = blue/purple
    limits    = c(0, 840),
    oob       = scales::squish,  # values > 840 shown as darkest colour
    na.value  = NA
  ) +
  theme_classic()
p_map_cwd


# ---- Map: gap-filled months -------------------------------------------------
# Separate palette family from the CWDmax map so the two are not confused. 
# Range 0-36: sum of gap-filled months per pixel across the
# three analysis years 
p_map_gap <- ggplot() +
  tidyterra::geom_spatraster(data = r_gap_filled_agg) +
  scico::scale_fill_scico(
    name      = "Gap-filled\nmonths [count]",
    palette   = "lajolla",
    direction = -1,              # low = light, high = dark red/brown
    limits    = c(0, 36),
    oob       = scales::squish,
    na.value  = NA
  ) +
  theme_classic()
p_map_gap


# ---- Sample for the distribution plot ---------------------------------------
# na.rm = FALSE plus drop_na() afterwards: na.rm = TRUE makes spatSample
# oversample repeatedly until enough valid pixels are found, which is
# prohibitively slow on a raster of this size.
set.seed(42)
df_sample <- terra::spatSample(
  r_conus,
  size   = 100000,
  method = "random",
  na.rm  = FALSE,
  as.df  = TRUE
) |>
  tibble::as_tibble() |>
  tidyr::drop_na()


# Gap-filled months shares
gap_shares <- data.frame(
  n_sampled      = nrow(df_sample),
  share_none     = round(100 * mean(df_sample$gap_filled_months == 0), 1),
  share_max_3    = round(100 * mean(df_sample$gap_filled_months <= 3), 1),
  share_above_12 = round(100 * mean(df_sample$gap_filled_months > 12), 1)
)

print(gap_shares)


# ---- Distribution: gap-filled months ----------------------------------------
# Discrete count (0-36) -> one bar per integer value, relative frequency.
p_bar_gap <- ggplot(df_sample, aes(x = gap_filled_months)) +
  geom_bar(
    aes(y = after_stat(count / sum(count))),
    fill = "grey30", width = 0.8
  ) +
  scale_x_continuous(breaks = seq(0, 36, 6)) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Gap-filled months (count)",
    y = "Share of pixels"
  ) +
  theme_classic()
p_bar_gap


# ---- Combined figure: gap-filled months -------------------------------------
# Map and distribution as one thesis figure. Titles are dropped here, the
# caption describes both panels.
p_gap_combined <- cowplot::plot_grid(
  p_map_gap + labs(title = NULL),
  p_bar_gap,
  ncol        = 1,
  labels      = "auto",
  align       = "v",
  axis        = "lr",          # panel edges stay flush despite the map legend
  rel_heights = c(1.6, 1)
)
p_gap_combined


# ---- Save -------------------------------------------------------------------
ggsave(here("fig", "map_cwd_max_conus.pdf"), p_map_cwd,
       width = 20, height = 14, units = "cm")
ggsave(here("fig", "map_cwd_max_conus.png"), p_map_cwd,
       width = 20, height = 14, units = "cm", dpi = 600)

ggsave(here("fig", "gap_filled_conus_combined.pdf"), p_gap_combined,
       width = 18, height = 16, units = "cm")
ggsave(here("fig", "gap_filled_conus_combined.png"), p_gap_combined,
       width = 18, height = 16, units = "cm", dpi = 600)