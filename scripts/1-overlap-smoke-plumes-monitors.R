source("script/functions.R")
source("script/packages.R")


# read in monitor locations
improve_loc <- read_csv("data/IMPROVE_CSN_sites.csv")
improve_loc_sf <- sfheaders::sf_point(improve_loc[,c("Longitude","Latitude")]) %>% mutate(Site = improve_loc$Site)

#read in plumes
smoke <- read_rds("data/hms/hms_smoke_plumes.rds")


#first do it for other plumes
ov <- list(); length(ov)<-length(smoke); names(ov)<-names(smoke)
for(i in 1:length(smoke)){

  if(nrow(smoke[[i]])>0){
    
  if("Density" %in% names(smoke[[i]]) == T){smoke[[i]]$density_missing = 0}
  if("Density" %in% names(smoke[[i]]) == F){smoke[[i]]$Density = NA; smoke[[i]]$density_missing = 1}
  
  smoke[[i]] <- smoke[[i]][st_is_valid(smoke[[i]]),]
  st_crs(improve_loc_sf) <- st_crs(smoke[[i]]) <- NA
  overlap <-st_join(improve_loc_sf, smoke[[i]], left = T) %>% 
    dplyr::filter(!is.na(Start) ) %>% 
    group_by(Site) %>% 
    summarise(low_count = sum(as.numeric(as.character(Density))==5, na.rm = T),
              med_count = sum(as.numeric(as.character(Density))==16, na.rm = T),
              high_count = sum(as.numeric(as.character(Density)) == 27, na.rm = T),
              density_missing = max(density_missing)) %>% 
    mutate(smoke_day = 1) #if there is overlap there is a smoke-day
  }else{overlap <- data.frame(Site =improve_loc_sf$Site, low_count = 0, med_count = 0, high_count = 0, density_missing = 1, smoke_day = 0, geometry  = NA )}
  
  
  
  
  ov[[i]] <- left_join(improve_loc, overlap) %>% 
    mutate(
      date = as.Date(names(smoke)[i], format = "%Y%m%d"),
      smoke_day = replace(smoke_day, is.na(smoke_day),0),
      low_count = replace(low_count, is.na(low_count),0),
      med_count = replace(med_count, is.na(med_count),0),
      high_count = replace(high_count, is.na(high_count),0),
      density_missing = replace(density_missing, is.na(density_missing),1)) %>% 
      dplyr::select(-geometry) %>% 
      dplyr::select(Network, Site, State, Latitude, Longitude, Elevation, date, smoke_day, low_count, med_count, high_count, density_missing) #density missing is wrong
  
  if(round(i/50)==i/50){print(paste("done with ",i," of ",length(smoke)))}
}

improve_smoke <- data.frame(data.table::rbindlist(ov)) 


smoke <- read_rds("data/hms/hms_smoke_plumes.rds")

density <- data.frame(date = as.Date(names(smoke), format = "%Y%m%d"), density_missing = NA)
for(i in 1:length(smoke)){
  
  if("Density" %in% names(smoke[[i]]) == F){density$density_missing[i] <- 1}
  if("Density" %in% names(smoke[[i]]) == T){density$density_missing[i] <- 1 - as.numeric(sum(!is.na(smoke[[i]]$Density) )>0) }
  
}


improve_smoke<-improve_smoke %>% dplyr::select(-density_missing)
improve_smoke <- left_join(improve_smoke, density)



write_rds(improve_smoke, file = "data/clean/improve_daily_hms_smoke.rds", compress = "gz")
