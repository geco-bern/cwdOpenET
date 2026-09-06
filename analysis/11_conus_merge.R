library(here)
library(terra)
library(ggplot2)
library(tidyterra)

# Read and merge to disk
files <- list.files(
  here("data", "conus_grid"),
  pattern    = "\\.tif$",
  full.names = TRUE
)

tiles <- lapply(files, terra::rast)
col   <- terra::sprc(tiles)

cwd <- terra::merge(
  col,
  filename  = here("data", "cwd_conus.tif"),
  overwrite = TRUE,
  datatype  = "FLT4S",
  gdal      = c("COMPRESS=DEFLATE", "TILED=YES", "BIGTIFF=YES")
)

# Sanity-Check: Resolution, CRS (EPSG:5070), bandnames
print(cwd)
names(cwd)
