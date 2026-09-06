# 09_h3_geology_edwards_landcover.R
# Extension of 08_h3_geology_edwards_plateau.R. checks whether
# the carbonate/clastic CWDmax contrast is hidden by the vegetation cover. 


# ---- Setup ------------------------------------------------------------------
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(readr)
library(here)
library(effectsize)

source(here("R", "add_boxplot_n.R"))
source(here("R", "lanid_filter.R"))
source(here("R", "save_fig.R"))

set.seed(42)

n_per_class <- 100000   # stratified sample per geology class (before veg split)


# Load CWD raster
r     <- terra::rast(here("data", "cwd_edwards_plateau.tif"))
r_cwd <- r[["cwd_max"]]


# Load, rasterize and classify geology 
sgmc <- sf::st_read(
  here("data", "geology_edwards_plateau.gpkg"),
  quiet = TRUE
) |>
  dplyr::filter(GENERALIZED_LITH %in% c(
    "Sedimentary, carbonate",
    "Sedimentary, clastic"
  )) |>
  dplyr::mutate(
    geo_class = dplyr::case_when(
      GENERALIZED_LITH == "Sedimentary, carbonate" ~ "Carbonate",
      GENERALIZED_LITH == "Sedimentary, clastic"   ~ "Clastic"
    ),
    geo_class = factor(geo_class, levels = c("Clastic", "Carbonate"))
  )

sgmc_vect           <- terra::vect(sgmc)
sgmc_vect$class_int <- as.integer(sgmc$geo_class)
geo_rast            <- terra::rasterize(sgmc_vect, r_cwd, field = "class_int")
names(geo_rast)     <- "geo_class"


# Load NLCD
stack_pre    <- terra::rast(here("data", "stack_edwards_plateau.tif"))
nlcd_aligned <- stack_pre[["nlcd"]]


# NLCD land cover shares, full Edwards Plateau extent
nlcd_labels <- c(
  `11` = "Open Water", `21` = "Developed, Open Space",
  `22` = "Developed, Low Intensity", `23` = "Developed, Medium Intensity",
  `24` = "Developed, High Intensity", `31` = "Barren Land",
  `41` = "Deciduous Forest", `42` = "Evergreen Forest", `43` = "Mixed Forest",
  `52` = "Shrub/Scrub", `71` = "Grassland/Herbaceous",
  `81` = "Pasture/Hay", `82` = "Cultivated Crops",
  `90` = "Woody Wetlands", `95` = "Emergent Herbaceous Wetlands"
)

nlcd_composition <- terra::freq(nlcd_aligned) |>
  dplyr::as_tibble() |>
  dplyr::transmute(code  = value,
                   class = nlcd_labels[as.character(value)],
                   n     = count,
                   share = round(n / sum(n), 4)) |>
  dplyr::arrange(dplyr::desc(share))

print(nlcd_composition, n = Inf)


# Three dominant natural/unmanaged classes from the table above.
# 1 = Grassland, 2 = Shrub, 3 = Forest
nlcd_reclass <- matrix(
  c(
    71, 1,          # Grassland/Herbaceous
    52, 2,          # Shrub/Scrub
    41, 3,          # Deciduous Forest
    42, 3,          # Evergreen Forest
    43, 3           # Mixed Forest
  ),
  ncol = 2, byrow = TRUE
)

veg_rast        <- terra::classify(nlcd_aligned, nlcd_reclass, others = NA)
names(veg_rast) <- "veg_class"


# Load LANID
lanid_aligned <- stack_pre[["lanid"]]


# Build combined class raster
# Mask geology to the three target vegetation classes and to pixels that
# LANID marks as non-irrigated.
combined <- terra::mask(geo_rast, veg_rast) # drop pixels outside the 3 veg classes
combined <- mask_irrigated(combined, lanid_aligned)
names(combined) <- "geo_class"


# ---- Stratified sample on geology class -------------------------------------
# n_per_class draws per geo_class stratum, not by area, so Clastic and
# Carbonate are comparably represented. veg_class is not itself stratified
sample_pts <- terra::spatSample(
  combined,
  size   = n_per_class,
  method = "stratified",
  na.rm  = TRUE,
  xy     = TRUE,
  as.df  = TRUE
)

ex_cwd <- terra::extract(r,        sample_pts[, c("x", "y")])
ex_veg <- terra::extract(veg_rast, sample_pts[, c("x", "y")])

veg_levels <- c("Grassland", "Shrub", "Forest")

df <- dplyr::tibble(
  geo_class         = sample_pts$geo_class,
  cwd_max           = ex_cwd$cwd_max,
  gap_filled_months = ex_cwd$gap_filled_months,
  veg_int           = ex_veg$veg_class
) |>
  dplyr::mutate(
    geo_class = factor(geo_class,
                       levels = c(1, 2),
                       labels = levels(sgmc$geo_class)),
    veg_class = factor(veg_int,
                       levels = seq_along(veg_levels),
                       labels = veg_levels)
  ) |>
  
  dplyr::filter(
    !is.na(cwd_max),
    !is.na(veg_class)
  )


# ---- Group statistics -------------------------------------------------------
stats_geo_veg <- df |>
  dplyr::group_by(veg_class, geo_class) |>
  dplyr::summarise(
    n       = dplyr::n(),
    median  = round(median(cwd_max), 1),
    .groups = "drop"
  )

# Substrate contrast within each vegetation class. Levels are
# c("Clastic", "Carbonate"), so a positive delta means Clastic holds the
# higher values.
contrasts_veg <- lapply(levels(df$veg_class), function(v) {
  d <- droplevels(dplyr::filter(df, veg_class == v))
  
  data.frame(
    veg_class   = v,
    diff_median = round(median(d$cwd_max[d$geo_class == "Clastic"]) -
                          median(d$cwd_max[d$geo_class == "Carbonate"]), 1),
    delta       = round(effectsize::cliffs_delta(cwd_max ~ geo_class, data = d,
                                                 ci = NULL)[[1]], 3)
  )
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(p_clastic_higher = round((delta + 1) / 2, 3))

# Type II ANOVA, since the vegetation classes are unequally distributed across the
# two substrates. The interaction term is the contrast varying by veg_class.
eta <- effectsize::eta_squared(
  car::Anova(lm(cwd_max ~ geo_class * veg_class, data = df), type = 2),
  partial = FALSE, ci = NULL
)

# Pooled across substrate. The sample is stratified by geo_class
stats_veg <- df |>
  dplyr::group_by(veg_class) |>
  dplyr::summarise(
    n       = dplyr::n(),
    median  = round(median(cwd_max), 1),
    .groups = "drop"
  )

# Vegetation composition within each substrate.
composition <- df |>
  dplyr::count(geo_class, veg_class) |>
  dplyr::group_by(geo_class) |>
  dplyr::mutate(share = round(n / sum(n), 3)) |>
  dplyr::ungroup()

print(as.data.frame(stats_geo_veg))
print(contrasts_veg)
print(as.data.frame(stats_veg))
print(as.data.frame(composition))
print(data.frame(parameter = eta$Parameter, eta2 = signif(eta$Eta2, 3)))


# Pooled comparison of clastic vs. carbonate but restricted to the three
# dominant vegetation classes
geo_cols <- c("Clastic" = "#0072B2", "Carbonate" = "#E69F00")

p_box_veg_subset <- ggplot(
  df,
  aes(x = geo_class, y = cwd_max, fill = geo_class)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2) +
  scale_fill_manual(values = geo_cols) +
  labs(
    x = NULL,
    y = expression(CWD[max]~"[mm]")) +
  theme_classic() +
  theme(legend.position = "none")

p_box_veg_subset <- p_box_veg_subset + add_boxplot_n(df, "geo_class", "cwd_max")

save_fig(p_box_veg_subset, "h3_edwards_boxplot_veg_subset", width = 12, height = 10)


# --- Faceted boxplot ---------------------------------------------------------
# Same Wong-palette colours as 09_h3_geology_edwards_plateau.R so that
# geo_class is always encoded identically across figures.

p <- ggplot(df, aes(x = veg_class, y = cwd_max, fill = geo_class)) +
  geom_boxplot(
    outlier.size  = 0.3,
    outlier.alpha = 0.2,
    position      = position_dodge(width = 0.8)
  ) +
  scale_fill_manual(values = geo_cols, name = NULL) +
  facet_wrap(~ veg_class, scales = "free_x", nrow = 1) +
  labs(
    x = NULL,
    y = expression(CWD[max]~"[mm]")
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    legend.position  = "bottom"
  )

# Per (veg_class x geo_class) n label.
df_n <- df |>
  dplyr::count(veg_class, geo_class)

p <- p +
  geom_text(
    data        = df_n,
    mapping     = aes(x = veg_class, y = max(df$cwd_max, na.rm = TRUE) * 1.05,
                      label = paste0("n = ", n), group = geo_class),
    position    = position_dodge(width = 0.8),
    inherit.aes = FALSE,
    size        = 2.6
  )

save_fig(p, "h3_edwards_landcover_boxplot", width = 16, height = 10)