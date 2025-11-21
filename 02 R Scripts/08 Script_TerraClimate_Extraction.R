
###############################################
### ENVIRONMENTAL DATA EXTRACTION - TerraClimate
### Averaged for each site (2018–2022)
###############################################

# ppt = precipitation (mm); pet = potential evapotranspiration (mm); tmax (monthly mean) = maximum temperature (ºC); tmin (monthly mean) = minimum temperature (ºC); vpd (monthly mean) = Vapor Pressure Deficit (kpd)

library(ncdf4)
library(raster)
library(readxl)
library(writexl)
library(dplyr)
library(purrr)
library(tibble)  # for tibble()

# Base folder containing TerraClimate rasters
tc_path <- "~/01 Masters_LA/07_environ_data_TerraClimate"

# Load site coordinates
sites <- read_excel("~/01 Masters_LA/04 Maps/moderadores.xlsx")

# Spatial extent to crop rasters (reduces memory usage)
roi_extent <- extent(-44.608007, -40.974161, -22.859079, -21.272083)

# Function:
# - loads raster stack
# - renames layers (month_year)
# - crops to extent
# - extracts raster values at site coordinates
# - calculates temporal mean (2018–2022)
extract_tc_mean <- function(variable, output_col) {
  
  # 1) Load TerraClimate rasters
  rst <- stack(list.files(
    tc_path,
    pattern = paste0("TerraClimate_", variable),
    full.names = TRUE
  ))
  
  # 2) Rename layers (1_2018, 2_2018, ... 12_2022)
  names(rst) <- apply(
    expand.grid(1:12, 2018:2022),
    1,
    FUN = paste,
    collapse = "_"
  )
  
  # 3) Crop raster
  rst_crop <- crop(rst, roi_extent)
  
  # 4) Extract climate values at each site
  raw_values <- extract(rst_crop, sites[, c("Longitude", "Latitude")])
  raw_df <- as.data.frame(raw_values)
  
  # 5) Temporal mean for each site (base R + tibble)
  mean_values <- rowMeans(raw_df, na.rm = TRUE)
  mean_df <- tibble(!!output_col := mean_values)
  
  return(mean_df)
}

# TerraClimate variables and output column names
tc_vars   <- c("ppt", "tmax", "tmin", "pet", "vpd")
tc_output <- c("mean_ppt", "mean_tmax", "mean_tmin", "mean_pet", "mean_vpd")

# Apply extraction function to each variable (bind columns)
tc_means <- purrr::map2_dfc(tc_vars, tc_output, extract_tc_mean)

# Final dataset: site info + climate means
env_data_final <- dplyr::bind_cols(sites, tc_means)

# Export final dataset
writexl::write_xlsx(
  env_data_final,
  "01 Datasets/02_processed_data/amb_TerraClimate_means_2018_2022.xlsx"
)


