# Emma Krasovich Southworth, emmars@stanford.edu
# originally created by Sam Heft-Neal, samhn@stanford.edu
# Last Updated: Jan 16, 2023
# Description: This function merges smoke plume data with the CSN + Improve Sites

# smoke_plumes_fp = file.path(wip_gdrive_fp, 'intermediate/hms_smoke_plumes.rds')
# loadd(c(spec_w_smoke_pm_df, CONUS_spec_sites_df), cache =drake_cache)

# function
merge_sites_w_smoke_plumes <- function(smoke_plumes_fp, spec_w_smoke_pm_df, CONUS_spec_sites_df) {


# CONVERT MONITOR DATAFRAME TO A SPATIAL OBJECT USING SF -----------------------
# turn monitor locations into an SF object:
  sites_sf <- CONUS_spec_sites_df %>% 
    st_as_sf(coords = c(x = 'long', y = 'lat'), crs = 4326)


# READ IN SMOKE DATA -----------------------------------------------------------
# read in smoke data:
smoke <- readRDS(smoke_plumes_fp) 

  # first do it for other plumes
  ov <- list(); length(ov)<-length(smoke); names(ov)<-names(smoke)
  
for(i in 1:length(smoke)){
  
  if(nrow(smoke[[i]])>0){
    
    if("Density" %in% names(smoke[[i]]) == T) {
      smoke[[i]]$density_missing = 0} # end if statement
    
    if("Density" %in% names(smoke[[i]]) == F) {
      smoke[[i]]$Density = NA; smoke[[i]]$density_missing = 1} # end if statement
    
    st_crs(smoke[[i]]) <- st_crs(sites_sf) <- NA # EK dont make this NA, keep it in a projection
    smoke[[i]] <- smoke[[i]][st_is_valid(smoke[[i]]),]
    
   # mapview(smoke[[i]], cex = 5) + mapview(overlap) + mapview(sites_sf, cex = 3)
    
    overlap <- st_join(sites_sf, smoke[[i]], left = T) %>% 
      dplyr::filter(!is.na(Start)) %>%
      group_by(site_id) %>%
      summarise(low_count = sum(as.numeric(as.character(Density))==5, na.rm = T),
                med_count = sum(as.numeric(as.character(Density))==16, na.rm = T),
                high_count = sum(as.numeric(as.character(Density))== 27, na.rm = T),
                density_missing = max(density_missing)) %>%
      mutate(smoke_day = 1) #if there is overlap there is a smoke-day
    
  }else{
    overlap <- data.frame(
      site_id = sites_sf$site_id, 
      low_count = 0,
      med_count = 0,
      high_count = 0,
      density_missing = 1, 
      smoke_day = 0, 
      geometry  = NA)
    } # end if else statement
  
  
  ov[[i]] <- left_join(spec_w_smoke_pm_df, overlap) %>% 
    mutate(
      Date = as.Date(names(smoke)[i], format = "%Y%m%d"),
      smoke_day = replace(smoke_day, is.na(smoke_day),0),
      low_count = replace(low_count, is.na(low_count),0),
      med_count = replace(med_count, is.na(med_count),0),
      high_count = replace(high_count, is.na(high_count),0),
      density_missing = replace(density_missing, is.na(density_missing),1)
      ) %>% 
    dplyr::select(-geometry) %>% 
    dplyr::select(Dataset,site_id, Date,
                  lat, long, epsg, smoke_day,  
                  low_count, med_count, high_count, density_missing) #density missing is wrong

  if(round(i/50)==i/50){print(paste("done with ",i," of ", length(smoke)))}
}

improve_smoke <- data.frame(data.table::rbindlist(ov)) 

# calculate density
density <- data.frame(Date = as.Date(names(smoke), format = "%Y%m%d"), 
                      density_missing = NA)
for(i in 1:length(smoke)){
  
  if("Density" %in% names(smoke[[i]]) == F){
    density$density_missing[i] <- 1}
  if("Density" %in% names(smoke[[i]]) == T){
    density$density_missing[i] <- 1 - as.numeric(sum(!is.na(smoke[[i]]$Density))>0) }
  
}

smoke_days_w_plumes <- improve_smoke %>% dplyr::select(-density_missing)
smoke_days_w_plumes <- left_join(smoke_days_w_plumes, density, by = 'Date')

return(smoke_days_w_plumes)

} # end function
