# Emma Krasovich Southworth, emmars@stanford.edu
# originally created by Sam Heft-Neal, samhn@stanford.edu
# Last Updated: Jan 16, 2023
# Description: This function merges smoke plume data with the CSN + Improve Sites

# smoke_plumes_fp = file.path(wip_gdrive_fp, 'intermediate/hms_smoke_plumes.rds')
# loadd(c(cleaned_spec_df, all_spec_sites), cache =drake_cache)

# function
merge_sites_w_smoke_plumes <- function(smoke_plumes_fp, cleaned_spec_df, all_spec_sites) {

# CONVERT MONITOR DATAFRAME TO A SPATIAL OBJECT USING SF -----------------------
# turn monitor locations into an SF object:
  aqs_sf <- cleaned_spec_df %>% 
    distinct(Dataset, SiteCode, Longitude, Latitude, State, epsg, Elevation) %>% 
    st_as_sf(coords = c(x = 'Longitude', y = 'Latitude'), crs = 4326)


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
    
    st_crs(smoke[[i]]) <- st_crs(aqs_sf) <- NA # EK dont make this NA, keep it in a projection
    smoke[[i]] <- smoke[[i]][st_is_valid(smoke[[i]]),]
    
   # mapview(smoke[[i]], cex = 5) + mapview(overlap) + mapview(aqs_sf, cex = 3)
    
    overlap <- st_join(aqs_sf, smoke[[i]], left = T) %>% 
      dplyr::filter(!is.na(Start)) %>%
      group_by(SiteCode) %>%
      summarise(low_count = sum(as.numeric(as.character(Density))==5, na.rm = T),
                med_count = sum(as.numeric(as.character(Density))==16, na.rm = T),
                high_count = sum(as.numeric(as.character(Density))== 27, na.rm = T),
                density_missing = max(density_missing)) %>%
      mutate(smoke_day = 1) #if there is overlap there is a smoke-day
    
  }else{
    overlap <- data.frame(
      SiteCode = aqs_sf$SiteCode, 
      low_count = 0,
      med_count = 0,
      high_count = 0,
      density_missing = 1, 
      smoke_day = 0, 
      geometry  = NA)
    } # end if else statement
  
  
  ov[[i]] <- left_join(all_spec_sites, overlap) %>% 
    mutate(
      Date = as.Date(names(smoke)[i], format = "%Y%m%d"),
      smoke_day = replace(smoke_day, is.na(smoke_day),0),
      low_count = replace(low_count, is.na(low_count),0),
      med_count = replace(med_count, is.na(med_count),0),
      high_count = replace(high_count, is.na(high_count),0),
      density_missing = replace(density_missing, is.na(density_missing),1)
      ) %>% 
    dplyr::select(-geometry) %>% 
    dplyr::select(Dataset,SiteCode, Date,
                  Latitude, Longitude, epsg, Elevation, smoke_day,  
                  low_count, med_count, high_count, density_missing) #density missing is wrong

  if(round(i/50)==i/50){print(paste("done with ",i," of ", length(smoke)))}
}

improve_smoke <- data.frame(data.table::rbindlist(ov)) 

# read in smoke data
smoke <- readRDS(smoke_plumes_fp) 

density <- data.frame(Date = as.Date(names(smoke), format = "%Y%m%d"), 
                      density_missing = NA)
for(i in 1:length(smoke)){
  
  if("Density" %in% names(smoke[[i]]) == F){
    density$density_missing[i] <- 1}
  if("Density" %in% names(smoke[[i]]) == T){
    density$density_missing[i] <- 1 - as.numeric(sum(!is.na(smoke[[i]]$Density))>0) }
  
}

improve_smoke <- improve_smoke %>% dplyr::select(-density_missing)
improve_smoke <- left_join(improve_smoke, density, by = 'Date')

return(improve_smoke)

} # end function
