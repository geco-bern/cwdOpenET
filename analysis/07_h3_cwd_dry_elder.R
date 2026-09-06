# 07_h3_cwd_dry_elder.R
# Validates whether the CWDmax contrast reported by Hahm et al. (2019)
# between Elder Creek (Coastal Belt) and Dry Creek (Central Belt) is
# recoverable in the 30 m CWDmax product.

# Extracts CWDmax, slope and northness for all pixels within the two
# catchment polygons and reports them next to a satellite map and a CWDmax
# map of the two sites.


# ---- Setup ------------------------------------------------------------------

library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(readr)
library(maptiles)
library(tidyterra)
library(ggspatial)
library(cowplot)
library(here)
library(effectsize)

source(here("R", "add_boxplot_n.R"))
source(here("R", "save_fig.R"))

# Panel geometry (cm). The maps sit in the top row, the boxplot in a shorter
# second row.
map_w_cm   <- 11
panel_h_cm <- 12
box_h_cm   <- 8

# No LANID filter here because the catchments do not feature irrigated pixels
r_cwd_full <- terra::rast(here("data", "cwd_northern_california.tif"))[["cwd_max"]]


# ---- Load hand-drawn catchment polygons -------------------------------------
# Outlines are in NAD83 CONUS Albers (parameter-identical to the raster's
# EPSG:5070), so this transform is effectively a no-op. It is done explicitly to
# guarantee CRS alignment before the extract.
poly_elder <- sf::st_read(here("data-raw", "elder_creek_polygon.shp"), quiet = TRUE) |>
  dplyr::mutate(site = "Elder Creek")
poly_dry <- sf::st_read(here("data-raw", "dry_creek_polygon.shp"), quiet = TRUE) |>
  dplyr::mutate(site = "Dry Creek")

polys <- dplyr::bind_rows(poly_elder, poly_dry) |>
  dplyr::select(site) |>
  sf::st_transform(crs = terra::crs(r_cwd_full))

# site is constant across each catchment
sf::st_agr(polys) <- "constant"


# ---- Terrain ----------------------------------------------------------------
# Load slope and northness
win <- terra::ext(terra::vect(polys)) + 1000

r_cwd <- terra::crop(r_cwd_full, win)

terrain_vars <- terra::rast(here("data", "stack_northern_california.tif"))[[c("slope", "northness")]] |>
  terra::crop(win)

r_all <- c(r_cwd, terrain_vars)


# ---- Extract all pixels falling within each polygon -------------------------
# A pixel is included if its centre falls inside the polygon
df <- terra::extract(r_all, terra::vect(polys)) |>
  dplyr::mutate(
    site = factor(polys$site[ID], levels = c("Elder Creek", "Dry Creek"))
  ) |>
  dplyr::filter(!is.na(cwd_max))

readr::write_csv(df, here("data", "h3_cwd_validation.csv"))


# ---- Group statistics -------------------------------------------------------
stats_sites <- df |>
  dplyr::group_by(site) |>
  dplyr::summarise(
    n         = dplyr::n(),
    median    = round(median(cwd_max), 1),
    iqr       = round(IQR(cwd_max), 1),
    slope     = round(median(slope, na.rm = TRUE), 1),
    northness = round(mean(northness, na.rm = TRUE), 2),
    .groups   = "drop"
  )

# Levels are c("Elder Creek", "Dry Creek"), so a positive delta means
# Elder Creek holds the higher values.
delta <- effectsize::cliffs_delta(cwd_max ~ site, data = df, ci = NULL)[[1]]

# Share of the total variance in CWDmax explained by catchment membership.
eta <- effectsize::eta_squared(
  aov(cwd_max ~ site, data = df),
  partial = FALSE, ci = NULL
)$Eta2

# Same comparison with slope and aspect held constant, using the gentle and
# south-facing classes of the aspect analysis.
stats_flat_south <- df |>
  dplyr::filter(slope >= 2, slope < 10, northness < -0.33) |>
  dplyr::group_by(site) |>
  dplyr::summarise(
    n       = dplyr::n(),
    median  = round(median(cwd_max), 1),
    .groups = "drop"
  )

print(as.data.frame(stats_sites))
print(as.data.frame(stats_flat_south))

cat("Median difference, all pixels:          ",
    stats_sites$median[1] - stats_sites$median[2], " mm\n",
    "Median difference, gentle south-facing: ",
    stats_flat_south$median[1] - stats_flat_south$median[2], " mm\n",
    "Cliff's delta:                          ", round(delta, 3), "\n",
    "P(Elder > Dry):                         ", round((delta + 1) / 2, 3), "\n",
    "Eta squared (site):                     ", round(eta, 3), "\n",
    sep = "")


# ---- Shared map elements ----------------------------------------------------
# Polygons are back-transformed to WGS84 only for display.
polys_wgs <- sf::st_transform(polys, crs = 4326)

# Label anchor points
labels_wgs <- polys |>
  sf::st_point_on_surface() |>
  sf::st_transform(crs = 4326)

# Bounding box of both catchments plus a margin (metres, EPSG:5070) so the
# polygons are not glued to the panel edge, then to WGS84 for the tiles.
margin_m <- 3000
box <- sf::st_bbox(polys)
box["xmin"] <- box["xmin"] - margin_m
box["xmax"] <- box["xmax"] + margin_m
box["ymin"] <- box["ymin"] - margin_m
box["ymax"] <- box["ymax"] + margin_m

box_wgs <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(box), 4326))
xmin <- box_wgs[["xmin"]]; xmax <- box_wgs[["xmax"]]
ymin <- box_wgs[["ymin"]]; ymax <- box_wgs[["ymax"]]

ext <- terra::ext(xmin, xmax, ymin, ymax)

# Vertical label nudge (map degrees), converted from ~1 cm on the saved
# figure using the y-extent and the panel height.
deg_per_cm     <- (ymax - ymin) / panel_h_cm
nudge_y_labels <- 0.012 + 0.5 * deg_per_cm


# ---- Backgrounds ------------------------------------------------------------
tiles <- get_tiles(
  x        = ext,
  provider = "Esri.WorldImagery",
  zoom     = 13,
  crop     = TRUE,
  project  = FALSE
)

# The tiles double as the target grid, so both maps sit on one identical
# geometry and the panels differ in nothing but the background values.
r_cwd_tiles <- terra::project(r_cwd_full, tiles, method = "bilinear")

# Fixed limits shared with the CONUS and zoom figures, so CWDmax reads the
# same across the paper.
scale_cwd <- scale_fill_viridis_c(
  name      = expression(CWD[max] ~ "[mm]"),
  direction = -1,
  limits    = c(0, 840),
  oob       = scales::squish,
  na.value  = NA
)


# ---- Panel builder ----------------------------------------------------------
# Everything except the background layer, its fill scale and the caption is
# shared, so both maps carry identical outlines, labels, annotations, breaks
# and theme.
build_site_panel <- function(layer_background, scale_fill = NULL,
                             caption = NULL, label_colour = "white") {
  ggplot() +
    layer_background +
    scale_fill +
    geom_sf(data = polys_wgs, fill = NA, colour = "red", linewidth = 0.4) +
    geom_sf_text(data = labels_wgs, aes(label = site),
                 colour = label_colour, size = 3, fontface = "bold",
                 nudge_y = nudge_y_labels,
                 fun.geometry = sf::st_zm) +
    coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), crs = 4326) +
    
    annotation_scale(
      location   = "bl",
      width_hint = 0.15,
      style      = "bar",
      bar_cols   = c("black", "white"),
      height     = grid::unit(0.12, "cm"),
      text_cex   = 0.8
    ) +
    
    annotation_north_arrow(
      location    = "br",
      which_north = "grid",
      style       = north_arrow_fancy_orienteering(
        fill     = c("white", "grey40"),
        line_col = "grey20",
        text_col = "white"
      ),
      height = unit(0.8, "cm"),
      width  = unit(0.8, "cm")
    ) +
    
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
    labs(x = "Longitude", y = "Latitude", caption = caption) +
    
    theme_classic(base_size = 13) +
    theme(
      axis.ticks        = element_line(colour = "grey30", linewidth = 0.4),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text         = element_text(size = 8),
      # hjust = 1 keeps each label's end anchored under its tick mark
      axis.text.x       = element_text(angle = 45, hjust = 1),
      axis.title        = element_text(size = 11),
      plot.title        = element_text(size = 14),
      plot.caption      = element_text(size = 6, colour = "grey40", hjust = 1),
      legend.title      = element_text(size = 9),
      legend.text       = element_text(size = 8),
      legend.key.width  = grid::unit(0.35, "cm"),
      # a colour bar comes out five times legend.key.height, so 0.6 gives 3 cm
      legend.key.height = grid::unit(0.6, "cm")
    )
}

# Panel a: satellite map
p_map <- build_site_panel(
  layer_background = geom_spatraster_rgb(data = tiles,
                                         maxcell = terra::ncell(tiles)),
  caption          = "\u00a9 Esri World Imagery"
)

# Panel b: same map on CWDmax 
p_cwd <- build_site_panel(
  layer_background = geom_spatraster(data = r_cwd_tiles, aes(fill = cwd_max),
                                     maxcell = terra::ncell(tiles)),
  scale_fill       = scale_cwd
)

# Panel c: boxplot 
p_box <- ggplot(df, aes(x = site, y = cwd_max)) +
  geom_boxplot(fill = "grey70", outlier.size = 0.3, outlier.alpha = 0.2) +
  labs(x = NULL, y = expression(CWD[max] ~ "[mm]")) +
  theme_classic() +
  theme(legend.position = "none")

p_box <- p_box + add_boxplot_n(df, "site", "cwd_max")

# Combine and save
# align/axis keep all panels the same size despite the CWDmax legend on panel
# b only, and put the boxplot in line with panel a. NULL leaves the
# bottom-right cell empty, the empty label suppresses its tag.
fig <- cowplot::plot_grid(
  p_map, p_cwd,
  p_box, NULL,
  ncol        = 2,
  labels      = c("a", "b", "c", ""),
  align       = "hv",
  axis        = "tblr",
  rel_heights = c(panel_h_cm, box_h_cm)
)

print(fig)

save_fig(fig, "h3_elder_dry_with_map",
         width  = 2 * map_w_cm,
         height = panel_h_cm + box_h_cm)