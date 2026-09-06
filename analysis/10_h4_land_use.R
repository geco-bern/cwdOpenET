# 10_h4_land_use.R
# Tests H4 on the High Plains region. Does CWDmax differ across
# land-use types?
# The classes are: Cropland irrigated, Cropland rainfed,
# Grassland and Shrubland.

# Grassland and Shrubland are restricted to rainfed pixels only, so that
# irrigation varies only within cropland. 


# Setup 
library(terra)
library(dplyr)
library(ggplot2)
library(here)
library(effectsize)

source(here("R", "add_boxplot_n.R"))
source(here("R", "lanid_filter.R"))
source(here("R", "save_fig.R"))

set.seed(42)


# Load data 
r <- terra::rast(here("data", "cwd_high_plains.tif"))

stack_pre <- terra::rast(here("data", "stack_high_plains.tif"))
nlcd      <- stack_pre[["nlcd"]]
lanid     <- stack_pre[["lanid"]]
slope     <- stack_pre[["slope"]]


# Build the four land-use classes 
nlcd_classes <- terra::classify(
  nlcd,
  rcl    = matrix(c(82, 2,
                    71, 3,
                    52, 4),
                  ncol = 2, byrow = TRUE),
  others = NA
)
class_rainfed <- mask_irrigated(nlcd_classes, lanid)

cropland_irrigated <- terra::ifel(nlcd == 82 & lanid == 1, 1, NA)

# spatSample() can only stratify one layer, so cover() is used to merge them into one. 
# taking the values of class_rainfed wherever cropland_irrigated is NA.
class_land_use <- terra::cover(cropland_irrigated, class_rainfed)
names(class_land_use) <- "class"


# Stratified sample 
# size is the number of samples per stratum, so all four classes end up
# comparably represented for the group comparison even though they have
# different areal shares
n_sample <- 100000

sample_pts <- terra::spatSample(
  class_land_use,
  size   = n_sample,
  method = "stratified",
  na.rm  = TRUE,
  xy     = TRUE,
  as.df  = TRUE
)

cwd_at_pts <- terra::extract(r, sample_pts[, c("x", "y")])

df_sample <- dplyr::tibble(
  x                 = sample_pts$x,
  y                 = sample_pts$y,
  class             = sample_pts$class,
  cwd_max           = cwd_at_pts$cwd_max,
  gap_filled_months = cwd_at_pts$gap_filled_months
) |>
  dplyr::mutate(
    class = factor(
      class,
      levels = c(1, 2, 3, 4),
      labels = c("Cropland irrigated", "Cropland rainfed",
                 "Grassland", "Shrubland")
    )
  ) |>
  dplyr::filter(!is.na(cwd_max))


# ---- Group statistics -------------------------------------------------------
stats_lu <- df_sample |>
  dplyr::group_by(class) |>
  dplyr::summarise(
    n       = dplyr::n(),
    median  = round(median(cwd_max), 1),
    .groups = "drop"
  )

# A positive delta means the first class of the pair holds the higher values.
lu_pairs <- list(
  c("Cropland irrigated", "Cropland rainfed"),
  c("Cropland rainfed",   "Grassland"),
  c("Cropland rainfed",   "Shrubland")
)

contrasts_lu <- lapply(lu_pairs, function(p) {
  d <- df_sample |>
    dplyr::filter(class %in% p) |>
    dplyr::mutate(class = factor(class, levels = p))
  
  data.frame(
    contrast    = paste(p[1], "-", p[2]),
    diff_median = round(median(d$cwd_max[d$class == p[1]]) -
                          median(d$cwd_max[d$class == p[2]]), 1),
    delta       = round(effectsize::cliffs_delta(cwd_max ~ class, data = d,
                                                 ci = NULL)[[1]], 3)
  )
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(p_first_higher = round((delta + 1) / 2, 3))

eta <- effectsize::eta_squared(
  aov(cwd_max ~ class, data = df_sample), ci = NULL
)$Eta2

print(as.data.frame(stats_lu))
print(contrasts_lu)
cat("Eta squared: ", signif(eta, 3), "\n", sep = "")

# Median slope across the whole region 
cat("Median slope, whole region: ",
    round(as.numeric(terra::global(slope, median, na.rm = TRUE)), 1),
    "\n", sep = "")


# Boxplot 
plot_land_use <- ggplot(
  data = df_sample,
  aes(x = class, y = cwd_max)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2, fill = "grey85") +
  labs(
    x = NULL,
    y = expression(CWD[max]~"[mm]")) +
  theme_classic() +
  theme(legend.position = "none")

# Add per-group n label.
plot_land_use <- plot_land_use + add_boxplot_n(df_sample, "class", "cwd_max")

save_fig(plot_land_use, "h4_land_use_irrigation", width = 14, height = 10)