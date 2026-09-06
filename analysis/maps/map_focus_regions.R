# map_focus_regions.R
# Maps of the four BA focus regions and the Northern California subset.
# Two combined 2 x 2 figures are written: one on Esri satellite imagery and
# one showing CWDmax.

# Extent and CRS are read from the exported cwd_max rasters.
# Maps are drawn in EPSG:5070, datum = 4326 keeps
# the axis labels in degrees.


# ---- Setup --------------------------------------------------------------
library(maptiles)
library(terra)
library(tidyterra)
library(ggplot2)
library(ggspatial)   # annotation_scale(), annotation_north_arrow()
library(cowplot)
library(sf)
library(here)
library(scales)
library(maps)

# Define output sizes
panel_w_cm <- 16
panel_h_cm <- 12

# Colour scale for the CWDmax version, shared by all panels and by the
# legend below. Same limits as Fig. 2 and Fig. 4.
scale_cwd <- scale_fill_viridis_c(
  name      = expression(CWD[max] ~ "[mm]"),
  direction = -1,
  limits    = c(0, 840),
  oob       = scales::squish,
  na.value  = NA
)

# --- Use function to build maps ----------------------------------------------
# Builds one focus-region map for every exported cwd_max
# raster. basemap = "tiles" draws Esri satellite imagery, basemap = "cwd"
# draws the CWDmax band of the raster itself. Everything else is identical
# between the two versions.
make_region_map <- function(raster_path, title = NULL, basemap = "tiles",
                            zoom = 10, show_title = TRUE, north_arrow = TRUE,
                            lon_break_width = 1, subset_bbox = NULL,
                            zoom_bbox = NULL, show_inset = TRUE) {
  
  basemap <- match.arg(basemap, c("tiles", "cwd"))
  
  # Extent and CRS from the exported raster
  r       <- terra::rast(raster_path)
  crs_map <- sf::st_crs(terra::crs(r))
  
  # unname() is required
  e    <- as.vector(terra::ext(r))
  xmin <- unname(e["xmin"]); xmax <- unname(e["xmax"])
  ymin <- unname(e["ymin"]); ymax <- unname(e["ymax"])
  
  if (basemap == "tiles") {
    
    # Tiles are requested for a slightly larger area than the panel. The
    # basemap is warped from Web Mercator into EPSG:5070, and without the
    # margin that warp can leave a thin empty strip along the panel edge.
    buffer_m <- 0.02 * max(xmax - xmin, ymax - ymin)
    
    box_tiles <- sf::st_as_sfc(
      sf::st_bbox(
        c(xmin = xmin - buffer_m, ymin = ymin - buffer_m,
          xmax = xmax + buffer_m, ymax = ymax + buffer_m),
        crs = crs_map
      )
    )
    
    # project = TRUE reprojects the tiles into the CRS of x (EPSG:5070), so
    # the data layer is the one thing never resampled.
    tiles <- get_tiles(
      x        = box_tiles,
      provider = "Esri.WorldImagery",
      zoom     = zoom,
      crop     = TRUE,
      project  = TRUE
    )
    
    tiles <- terra::crop(tiles, terra::ext(xmin, xmax, ymin, ymax))
    
    main_map <- ggplot() +
      geom_spatraster_rgb(data = tiles)
    
  } else {
    
    # maxcell caps the number of cells drawn; lower it if plotting is slow.
    main_map <- ggplot() +
      geom_spatraster(data = r[[1]], maxcell = 5e5) +
      scale_cwd
    
  }
  
  # Subset outline: drawn after the basemap so it sits on the map.
  if (!is.null(subset_bbox)) {
    main_map <- main_map +
      annotate(
        "rect",
        xmin = subset_bbox[["xmin"]], xmax = subset_bbox[["xmax"]],
        ymin = subset_bbox[["ymin"]], ymax = subset_bbox[["ymax"]],
        fill = NA, colour = "red", linewidth = 0.9
      )
  }
  
  # Outline of the zoom window shown in Fig. 4, same style as the subset.
  if (!is.null(zoom_bbox)) {
    main_map <- main_map +
      annotate(
        "rect",
        xmin = zoom_bbox[["xmin"]], xmax = zoom_bbox[["xmax"]],
        ymin = zoom_bbox[["ymin"]], ymax = zoom_bbox[["ymax"]],
        fill = NA, colour = "red", linewidth = 0.9
      )
  }
  
  main_map <- main_map +
    # datum = 4326 keeps the axis labels in degrees. expand is left at its default (TRUE) so
    # ggplot pads a visible gap between the axis lines and the map
    coord_sf(
      crs   = crs_map,
      datum = sf::st_crs(4326),
      xlim  = c(xmin, xmax),
      ylim  = c(ymin, ymax)
    ) +
    
    annotation_scale(
      location   = "bl",
      width_hint = 0.15,
      style      = "bar",
      bar_cols   = c("black", "white"),
      height     = grid::unit(0.12, "cm"),
      text_cex   = 0.8
    ) +
    
    scale_x_continuous(breaks = scales::breaks_width(lon_break_width)) +
    labs(
      x = "Longitude",
      y = "Latitude"
    ) +
    
    theme_classic(base_size = 13) +
    theme(
      axis.ticks        = element_line(colour = "grey30", linewidth = 0.4),
      axis.ticks.length = unit(0.15, "cm"),
      axis.text         = element_text(size = 8),
      # Rotated longitude labels take less horizontal space per label.
      # hjust = 1 keeps the label's end anchored under its tick mark.
      axis.text.x       = element_text(angle = 45, hjust = 1),
      axis.title        = element_text(size = 11),
      plot.title        = element_text(size = 14),
      # One shared legend is added to the combined figure instead.
      legend.position   = "none"
    )
  
  if (show_title) {
    main_map <- main_map + labs(title = title)
  }
  
  if (north_arrow) {
    main_map <- main_map +
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
      )
  }
  
  # World map for inset. Skipped entirely when show_inset = FALSE, so
  # maps::map() isn't called for maps that don't need it.
  if (show_inset) {
    world <- sf::st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))
    
    # The study box is drawn in the WGS84 inset, where the EPSG:5070
    # rectangle is a tilted quadrilateral. st_segmentize() adds vertices
    # along the edges before transforming, otherwise only the corners move
    # and the curved edges are drawn as straight lines.
    study_box <- sf::st_as_sfc(
      sf::st_bbox(c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
                  crs = crs_map)
    ) |>
      sf::st_segmentize(dfMaxLength = 5000) |>
      sf::st_transform(4326)
    
    inset_map <- ggplot() +
      geom_sf(data = world, fill = "grey30", colour = "grey50",
              linewidth = 0.1) +
      geom_sf(data = study_box, fill = "red", colour = "red", alpha = 0.4) +
      coord_sf(xlim = c(-130, -60), ylim = c(20, 55), expand = FALSE) +
      theme_void(base_size = 7) +
      theme(
        panel.background = element_rect(fill = "lightblue", colour = NA),
        panel.border     = element_rect(fill = NA, colour = "black",
                                        linewidth = 0.8)
      )
    
    final_map <- ggdraw(main_map) +
      draw_plot(inset_map, x = 0.71, y = 0.75, width = 0.20, height = 0.20)
  } else {
    final_map <- ggdraw(main_map)
  }
  
  final_map
}


# ---- Generate maps ----------------------------------------------------------
# Zoom windows of Fig. 4, copied from map_zooms.R. The Edwards Plateau has
# no zoom window.
bbox_zoom_norcal      <- c(xmin = -2196800.0, xmax = -2060790.0,
                           ymin =  2143053.0, ymax =  2270823.0)
bbox_zoom_highplains  <- c(xmin =  -498940.6, xmax =  -376923.9,
                           ymin =  1597756.5, ymax =  1708584.7)
bbox_zoom_appalachia  <- c(xmin =  1134811.7, xmax =  1265610.5,
                           ymin =  1653744.2, ymax =  1774527.1)

# Northern California subset used for H2.
bbox_subset_norcal <- c(xmin = -2336730, xmax = -2229380,
                        ymin =  2138000, ymax =  2366160)

# Panel order for both figures. The satellite version carries the subset
# outline, the CWDmax version the zoom windows, so a red rectangle has one
# meaning per figure.
make_all_regions <- function(basemap) {
  tiles <- basemap == "tiles"
  cwd   <- basemap == "cwd"
  
  list(
    make_region_map(
      raster_path     = here("data", "cwd_northern_california.tif"),
      title           = "Northern California",
      basemap         = basemap,
      lon_break_width = 0.5,
      subset_bbox     = if (tiles) bbox_subset_norcal,
      zoom_bbox       = if (cwd)   bbox_zoom_norcal
    ),
    make_region_map(
      raster_path = here("data", "cwd_edwards_plateau.tif"),
      title       = "Edwards Plateau",
      basemap     = basemap
    ),
    make_region_map(
      raster_path = here("data", "cwd_high_plains.tif"),
      title       = "High Plains",
      basemap     = basemap,
      zoom_bbox   = if (cwd) bbox_zoom_highplains
    ),
    make_region_map(
      raster_path = here("data", "cwd_appalachia.tif"),
      title       = "Appalachians",
      basemap     = basemap,
      zoom_bbox   = if (cwd) bbox_zoom_appalachia
    )
  )
}

maps_tiles <- make_all_regions("tiles")
maps_cwd   <- make_all_regions("cwd")


# ---- Combined figure on satellite imagery -----------------------------------
fig_regions <- plot_grid(
  plotlist = maps_tiles,
  labels   = c("a", "b", "c", "d"),
  ncol     = 2
)

fig_regions <- ggdraw(fig_regions) +
  draw_label(
    "\u00a9 Esri World Imagery",
    x = 0.99, y = 0.01, hjust = 1, vjust = 0,
    size = 6, colour = "grey40"
  )


# ---- Combined figure showing CWDmax -----------------------------------------
fig_regions_cwd <- plot_grid(
  plotlist = maps_cwd,
  labels   = c("a", "b", "c", "d"),
  ncol     = 2
)

# The legend is built from a throwaway plot so that the four panels keep
# exactly the size they have in the satellite version.
legend_cwd <- get_plot_component(
  ggplot(data.frame(x = 1, y = 1, z = c(0, 840)), aes(x, y, fill = z)) +
    geom_raster() +
    scale_cwd +
    theme(legend.position  = "bottom",
          legend.key.width = unit(1.4, "cm")),
  "guide-box-bottom"
)

fig_regions_cwd <- plot_grid(fig_regions_cwd, legend_cwd,
                             ncol = 1, rel_heights = c(1, 0.05))


# ---- Display ----------------------------------------------------------------
fig_regions
fig_regions_cwd


# ---- Save -------------------------------------------------------------------
ggsave(here("fig", "map_focus_regions_combined.pdf"), plot = fig_regions,
       width = 2 * panel_w_cm, height = 2 * panel_h_cm, units = "cm",
       dpi = 600)

ggsave(here("fig", "map_focus_regions_combined.png"), plot = fig_regions,
       width = 2 * panel_w_cm, height = 2 * panel_h_cm, units = "cm",
       dpi = 600)

ggsave(here("fig", "map_focus_regions_cwd_combined.pdf"), plot = fig_regions_cwd,
       width = 2 * panel_w_cm, height = 2 * panel_h_cm + 1.5, units = "cm",
       dpi = 600)

ggsave(here("fig", "map_focus_regions_cwd_combined.png"), plot = fig_regions_cwd,
       width = 2 * panel_w_cm, height = 2 * panel_h_cm + 1.5, units = "cm",
       dpi = 600)