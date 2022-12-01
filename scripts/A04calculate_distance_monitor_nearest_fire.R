# Description: function to calculate the distance from fire to nearest monitoring station

# fires_fp = file.path(wip_gdrive_fp, 'intermediate/hms_fires.RDS')
# loadd(c(improve_smoke_dens_df, all_spec_sites), cache = drake_cache)


calculate_distance_monitor_nearest_fire <- function(fires_fp, improve_smoke_dens_df,
                                                    all_spec_sites) {

# turn monitor locations into sf object
aqs_sf <- all_spec_sites %>% 
    st_as_sf(crs = 4326, coords = c('Longitude', 'Latitude')) 
  

# read in fire data
fire_df <- read_rds(fires_fp) 

i  <- 1
firedist <- list()

for(i in 1:length(fire_df)){
  if(nrow(fire_df[[i]])>0){
    # set crs
    fire_df[[i]] <- fire_df[[i]][st_is_valid(fire_df[[i]]),]
    st_crs(fire_df[[i]]) <- st_crs(aqs_sf) <- NA
    
    fire_df[[i]] <- fire_df[[i]][unlist(lapply(fire_df[[i]]$geometry, function(x){is.nan(x)[1]}))==F,]
  
  test <- st_nn(x = aqs_sf, y = fire_df[[i]],k = 1, sparse = T, returnDist = T)[[2]] %>% unlist()


    firedist[[i]] <- data.frame(all_spec_sites, 
                            km2fire = round(test*111), 
                            Date = as.Date(names(fire_df)[i], format = "%Y%m%d"))
  }else{
    firedist[[i]] <- data.frame(all_spec_sites, 
                                km2fire = 1e6, 
                                Date = as.Date(names(fire_df)[i], format = "%Y%m%d"))
              } # end ifelse
} # end for loop

firedist <- data.frame(rbindlist(firedist))

# combine with improve smoke density data
improve_smoke_dens_fire_df <- improve_smoke_dens_df %>% 
  left_join(firedist %>% 
              dplyr::select(-c(Longitude, Latitude, epsg, Elevation)), 
            by = c('Dataset', 'SiteCode', 'Date')) %>% 
  distinct()

return(improve_smoke_dens_fire_df)

} # end function


