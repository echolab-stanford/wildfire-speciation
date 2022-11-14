# ------------------------------------------------------------------------------
# Merge "plumes_speciation_at_sites.csv" with wildfire PM2.5
# Written by: Ayako Kawano
# Last edited: Nov 2022
# ------------------------------------------------------------------------------

#### Settings ####
data_path <- "G:/Shared drives/echolab data/wildfire_speciation/intermediate/"

#### Load libraries ####
library(tools)
library(lubridate)
library(sf)
library(tigris)
library(rNOMADS)
library(FNN)
library(parallel)
library(dplyr)
library(tidyr)
library(purrr)
library(Metrics)
library(fixest)
library(ggplot2)
library(ggpmisc)
library(ncdf4)
library(sf)
library(raster)
library(RColorBrewer)
library(maptools)
library(spdplyr)
library(stringr)
library(rnaturalearth)
library(tidyverse)
library(RColorBrewer)
library(rgeos)
library(s2)
library(anytime)

#### Load data ####
# speciation data
df <- read.csv(paste(data_path, "plumes_speciation_at_sites.csv", sep = ''))

# wildfire pm data
pm <- readRDS("C:/Users/ayako/BurkeLab Dropbox/Projects/daily-10km-smokePM/final/10km_grid/smokePM2pt5_predictions_daily_10km_20060101-20201231.rds")
# Load 10 km grid
grid_10km <- read_sf("C:/Users/ayako/BurkeLab Dropbox/Projects/daily-10km-smokePM//final/10km_grid/10km_grid_wgs84/10km_grid_wgs84.shp")


#### Match monitor locations in speciation data with grid 10km ID ####
# make points spatial data for monitors in df
points <- st_as_sf(df,coords = c('long',"lat")) %>% dplyr::select(site_code, geometry)
# drop duplicates
points <- points %>% distinct() 
st_crs(points)= 4326 # same CRS as grid 10km
points$grid_id_10km <- apply(st_intersects(grid_10km, points, sparse = FALSE), 2, 
                                function(col) {grid_10km[which(col), ]$ID}) 

# Add grid_id_10km to the speciation df
speciation <- left_join(df, points, by = "site_code")

# Merge files to get pm prediction and add to speciation df
speciation$date <- anydate(speciation$date)
speciation$grid_id_10km <- as.integer(speciation$grid_id_10km)
speciation <- left_join(speciation, pm, by = c("grid_id_10km", "date"))

write_rds(speciation, paste(data_path, "smokePM_plumes_speciation_at_sites.rds", sep = ''))
write_csv(speciation, paste(data_path, "smokePM_plumes_speciation_at_sites.csv", sep = ''))
