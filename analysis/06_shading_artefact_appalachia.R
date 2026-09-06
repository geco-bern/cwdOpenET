# 06_shading_artefact_appalachia.R
# This script works in the same way as the 05_shading_artefact_northern_california.R one but
# for the appalachians focus region. 


# ---- Setup ------------------------------------------------------------------
library(terra)
library(tidyr)
library(dplyr)
library(ggplot2)
library(here)

source(here("R", "annotate_lm.R"))
source(here("R", "add_boxplot_n.R"))
source(here("R", "lanid_filter.R"))
source(here("R", "save_fig.R"))

set.seed(42)
n_sample <- 1000000

aspect_colours <- c(
  "south-facing" = "#d73027",
  "east/west"    = "#fee08b",
  "north-facing" = "#4575b4"
)

# Load data
cwd       <- terra::rast(here("data", "cwd_appalachia.tif"))[["cwd_max"]]

stack_pre <- terra::rast(here("data", "stack_appalachia.tif"))
elevation <- stack_pre[["elevation"]]
twi       <- stack_pre[["twi"]]
northness <- stack_pre[["northness"]]
slope     <- stack_pre[["slope"]]
nlcd      <- stack_pre[["nlcd"]]    
lanid     <- stack_pre[["lanid"]]


# ---- Sample -----------------------------------------------------------------
stack <- c(cwd, elevation, twi, northness, slope)
names(stack) <- c("cwd_max", "elevation", "twi", "northness", "slope")

df_raw <- terra::spatSample(
  stack,
  size   = n_sample * 2,
  method = "regular",
  na.rm  = FALSE,
  xy     = TRUE,
  as.df  = TRUE
) |>
  tibble::as_tibble()

# LANID value per sampled pixel, used by the drop_irrigated() call below.
df_raw$lanid <- terra::extract(lanid, as.matrix(df_raw[, c("x", "y")]))$lanid

# NLCD land cover, extracted here for the forest filter used in Plots 4 and 5.
df_raw$nlcd  <- terra::extract(nlcd, as.matrix(df_raw[, c("x", "y")]))$nlcd

# Derive classification variables used in the plots below:
# aspect_class buckets northness into south-/east-west-/north-facing
# terciles (breaks at ±0.33). slope_class groups slope into gentle/moderate/steep bins
# elev_band splits elevation into deciles to relate slope to elevation.

# Minimum slope for a meaningful aspect. On near-flat terrain aspect is
# undefined and the algorithm assigns a random value, so northness
# carries no information there.
slope_min <- 2   # degrees

# Slope class labels are defined once here and reused by every scale that
# refers to them
slope_labels <- c(
  gentle   = paste0("gentle (", slope_min, " to 10 deg)"),
  moderate = "moderate (10 to 25 deg)",
  steep    = "steep (>25 deg)"
)

df <- df_raw |>
  drop_irrigated() |>
  tidyr::drop_na(cwd_max, elevation, twi, northness, slope) |>
  dplyr::filter(slope >= slope_min) |>
  mutate(
    aspect_class = cut(
      northness,
      breaks         = c(-1, -0.33, 0.33, 1),
      labels         = c("south-facing", "east/west", "north-facing"),
      include.lowest = TRUE
    ),
    slope_class = cut(
      slope,
      breaks         = c(slope_min, 10, 25, 90),
      labels         = unname(slope_labels),
      include.lowest = TRUE
    ),
    elev_band = cut(
      elevation,
      breaks         = quantile(elevation, probs = seq(0, 1, 0.1), na.rm = TRUE),
      labels         = paste0("Band ", 1:10),
      include.lowest = TRUE
    )
  )

# Sample size per aspect class, used to annotate Plots 2 and 3 below.
label_n_aspect_appalachia <- df |>
  dplyr::count(aspect_class) |>
  dplyr::mutate(line = paste0(aspect_class, ": n = ", format(n, big.mark = "'"))) |>
  dplyr::pull(line) |>
  paste(collapse = "\n")


# ---- Plot 1: CWDmax vs northness --------------------------------------------
mod_northness_appalachia <- lm(cwd_max ~ northness, data = df)

plot_northness_appalachia <- ggplot(
  data = df,
  aes(x = northness, y = cwd_max)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "n (log10)") +
  geom_smooth(method = "lm", formula = y ~ x,
              se = FALSE, colour = "red", linewidth = 0.9) +
  scale_x_continuous(
    breaks       = c(-1, -0.5, 0, 0.5, 1),
    labels       = c("-1\nsouth", "-0.5", "0\neast/west", "0.5", "1\nnorth"),
    minor_breaks = NULL
  ) +
  labs(
    x = "Northness [cos(aspect)]",
    y = expression(CWD[max]~"[mm]")
  ) +
  theme_classic()

# Add R²/n/slope.
plot_northness_appalachia <- annotate_lm(
  plot_northness_appalachia, mod_northness_appalachia,
  slope_unit = "mm/unit northness"
)

save_fig(plot_northness_appalachia, "h2_appalachia_cwd_by_northness")


# ---- Plot 2: TWI vs CWDmax by aspect ----------------------------------------
# One lm per aspect class, matching the three geom_smooth() lines drawn below,
# so slope and R² can be reported per class.
fits_aspect_appalachia <- df |>
  dplyr::group_by(aspect_class) |>
  dplyr::group_modify(function(.x, .y) {
    m <- lm(cwd_max ~ twi, data = .x)
    tibble::tibble(
      n         = nrow(.x),
      slope     = round(coef(m)[["twi"]], 2),
      r_squared = round(summary(m)$r.squared, 4)
    )
  }) |>
  dplyr::ungroup()

print(as.data.frame(fits_aspect_appalachia))

plot_aspect_appalachia <- ggplot(
  data = df, 
  aes(x = twi, y = cwd_max, colour = aspect_class)) +
  geom_point(alpha = 0.04, size = 0.4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  scale_colour_manual(values = aspect_colours, name = "Aspect") +
  labs(
    x = "TWI",
    y = expression(CWD[max]~"[mm]")) +
  theme_classic() +
  theme(legend.position = "top")

# Add per-aspect-group n label 
plot_aspect_appalachia <- plot_aspect_appalachia +
  annotate(
    "label",
    x = Inf, y = Inf, hjust = 1.05, vjust = 1.1,
    label = label_n_aspect_appalachia,
    size = 3.2, label.size = 0.3,
    label.padding = unit(0.4, "lines"),
    fill = scales::alpha("white", 0.9)
  )

save_fig(plot_aspect_appalachia, "h2_appalachia_cwd_aspect")

# ---- Plot 3: slope x aspect ------------------------------------------------
# Fixed y-axis limits, shared with Plot 3 in
# 05_shading_artefact_northern_california.R for cross-region comparability.
shared_y_limits_slope_aspect <- c(137, 826)

# Mean CWDmax per aspect x slope class. These nine values are the points drawn
# below 
means_slope_aspect <- df |>
  group_by(aspect_class, slope_class) |>
  summarise(
    n        = dplyr::n(),
    mean_cwd = round(mean(cwd_max, na.rm = TRUE), 1),
    .groups  = "drop"
  )

print(as.data.frame(means_slope_aspect))

plot_slope_aspect <- means_slope_aspect |>
  ggplot(
    aes(x = slope_class, y = mean_cwd, colour = aspect_class, group = aspect_class)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_colour_manual(values = aspect_colours, name = "Aspect") +
  scale_x_discrete(limits = rev) +
  coord_cartesian(ylim = shared_y_limits_slope_aspect) +
  labs(
    x       = "Slope class",
    y       = expression(CWD[max]~"[mm]")) +
  theme_classic() +
  theme(legend.position = "top")

# Per-aspect-group n label 
plot_slope_aspect <- plot_slope_aspect +
  annotate(
    "label",
    x = Inf, y = Inf, hjust = 1.05, vjust = 1.1,
    label = label_n_aspect_appalachia,
    size = 3.2, label.size = 0.3,
    label.padding = unit(0.4, "lines"),
    fill = scales::alpha("white", 0.9)
  )

save_fig(plot_slope_aspect, "h2_appalachia_slope_aspect")


# ---- Plot 4: TWI vs CWDmax, south-facing forest only ------------------------
df_forest <- df |>
  dplyr::filter(nlcd %in% c(41, 42, 43)) |>   # deciduous/evergreen/mixed forest
  dplyr::filter(aspect_class == "south-facing")

# Plot
mod_twi_forest_south_appalachia <- lm(cwd_max ~ twi, data = df_forest)

plot_twi_forest_south <- ggplot(
  data = df_forest,
  aes(x = twi, y = cwd_max)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "n (log10)") +
  geom_smooth(method = "lm", formula = y ~ x,
              se = FALSE, colour = "red", linewidth = 0.9) +
  geom_smooth(method = "gam", formula = y ~ s(x), colour = "black",
              linewidth = 0.9, linetype = "dashed") +
  labs(
    x = "TWI",
    y = expression(CWD[max]~"[mm]")) +
  theme_classic()

# Add R²/n/slope label
plot_twi_forest_south <- annotate_lm(
  plot_twi_forest_south, mod_twi_forest_south_appalachia,
  slope_unit = "mm/TWI unit"
)

save_fig(plot_twi_forest_south, "h2_appalachia_twi_forest_south")


# ---- Plot 5: TWI vs CWDmax, forest only, all aspects pooled -----------------
df_forest_pooled <- df |>
  dplyr::filter(nlcd %in% c(41, 42, 43))   # deciduous/evergreen/mixed forest

# Plot
mod_twi_forest_pooled_appalachia <- lm(cwd_max ~ twi, data = df_forest_pooled)

plot_twi_forest_pooled <- ggplot(
  data = df_forest_pooled,
  aes(x = twi, y = cwd_max)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "n (log10)") +
  geom_smooth(method = "lm", formula = y ~ x,
              se = FALSE, colour = "red", linewidth = 0.9) +
  geom_smooth(method = "gam", formula = y ~ s(x), colour = "black",
              linewidth = 0.9, linetype = "dashed") +
  labs(
    x = "TWI",
    y = expression(CWD[max]~"[mm]")) +
  theme_classic()

plot_twi_forest_pooled <- annotate_lm(
  plot_twi_forest_pooled, mod_twi_forest_pooled_appalachia,
  slope_unit = "mm/TWI unit"
)

save_fig(plot_twi_forest_pooled, "h2_appalachia_twi_forest_pooled")

# Median slope for Appalachian region
cat("Median slope, region: ",
    round(as.numeric(terra::global(slope, median, na.rm = TRUE)), 1),
    "\n", sep = "")
