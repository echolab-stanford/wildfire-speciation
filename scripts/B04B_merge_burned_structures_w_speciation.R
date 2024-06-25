# Jeff Wen, Emma Krasovich Southworth
# Last Updated: April 2, 2024 by EKS
# Description: estimating relationship of chemical species concentration for fires that burned structures
# and smoke

# inputs:
# loadd(c(globfire_structure_joined_df,CONUS_spec_sites_df, clean_PMspec_df, clean_PMspec_sites_df), cache = drake_cache)
# grid_fp = file.path(data_fp, 'intermediate/10km_grid_wgs84.shp')

# FUNCTION
merge_burned_structures_w_speciation <- function(globfire_structure_joined_df,clean_PMspec_sites_df, 
                                                 clean_PMspec_df, grid_fp) {
  # parellelize
  future::plan(multisession) 
  
# Monitor, GlobFire merge ----
## load grid cell data
grid_transform <- st_read(grid_fp) %>%
  st_transform(crs=4326)

## monitor location lat lon merge to get gridcell id
monitor_location_df <- clean_PMspec_sites_df %>%
  st_as_sf(coords=c("long","lat"), crs=4326, agr="constant") %>%
  # st_transform("epsg:5070") %>%
  st_join(grid_transform)

## merge monitor location/gridcell with readings data
monitor_reading_df <- clean_PMspec_df %>%
  dplyr::select(Dataset, state_name, year, month, Date, site_id, nonsmokePM_MF,
                MF_adj, smokePM, AL:ZR, lat, long) %>%
  left_join(monitor_location_df %>%
              st_drop_geometry() %>%
              dplyr::select(site_id, ID, COORDX, COORDY),
            by=c("site_id"))

# get all smokepm files at 10km
fire_files_list <- list.files(
  file.path(data_fp, 'intermediate/fire_smokepm/30km'),
  full.names = TRUE)

# current_fp <- fire_files_list[1]

gridcell_smokepm_dt <- purrr::map_df(fire_files_list, function(current_fp){

  temp_fire_smokepm_dt <- read_fst(current_fp,
                                   as.data.table = T) %>% 
    mutate(smokePM_pred_coded = ifelse(is.na(smokePM_pred_coded), 0, smokePM_pred_coded))

  # filter smokePM to the monitor location gridcells and where there is smokePM
  # we want to compare smoke speciation for fires with specific characteristics vs without
  # so to be sure we filter for fires matched to specific fire_ids (there are cells with smokepm not matched to fire_ids...)
  temp_fire_smokepm_filt_dt <- temp_fire_smokepm_dt[
    ID %in% unique(monitor_reading_df$ID) & smokePM_pred_coded>=0 & !is.na(fire_id),
  ][,`:=`(month=month(date),
          year=year(date))
  ]
  return(temp_fire_smokepm_filt_dt)
}) %>%
  bind_rows()

# Globfire, structures merge ----
## ending with monitor-contributedpm/#structures-day
grid_smoke_structure_dt <- gridcell_smokepm_dt %>%
  left_join(globfire_structure_joined_df %>%
              mutate(Id = as.character(Id),
                     fire_duration=as.integer(FDate-IDate)) %>%
              mutate(fire_duration = ifelse(fire_duration == 0, 1, fire_duration), # for fires that were only 1 day, calculation gives 0, but should be 1
                     daily_structures_destroyed=structures_destroyed/fire_duration) %>%
              dplyr::select(Id, year, IDate, FDate, fire_duration, incid_name_lower, coverage_perc,
                            structures_destroyed, daily_structures_destroyed),
            by=c("fire_id"="Id", "year"="year")) %>%
  mutate(structures_destroyed=ifelse(is.na(structures_destroyed), 0, structures_destroyed),
         contrib_structures_destroyed=share*structures_destroyed,
         daily_structures_destroyed=ifelse(is.na(daily_structures_destroyed), 0, daily_structures_destroyed),
         contrib_daily_structures_destroyed=share*daily_structures_destroyed) %>%
  group_by(ID, month, year, date) %>%
  dplyr::summarize(avg_cum_traj_distance=mean(avg_cum_traj_distance, na.rm=T),
            smokePM_pred=mean(smokePM_pred, na.rm=T),
            smokePM_pred_coded=mean(smokePM_pred_coded, na.rm=T),
            contrib_smokePM=sum(contrib_smokePM, na.rm=T),
            structures_destroyed=mean(structures_destroyed, na.rm=T),
            contrib_structures_destroyed=sum(contrib_structures_destroyed, na.rm=T),
            daily_structures_destroyed=sum(daily_structures_destroyed, na.rm=T),
            contrib_daily_structures_destroyed=sum(contrib_daily_structures_destroyed, na.rm=T)) %>% 
  ungroup()

## aggregate by monitor-day to calc weighted avg. of structures burned, weighted avg. of soil characteristics
merged_smoke_spec_burned_structures_df <- monitor_reading_df %>% 
  left_join(grid_smoke_structure_dt %>%
              dplyr::select(ID, month, year, Date ='date', contrib_smokePM, contrib_daily_structures_destroyed),
            by=c("ID","month","year","Date")) #%>% 
  # filter(!is.na(contrib_daily_structures_destroyed)) %>% 
  # dplyr::select(ID, Dataset:site_id, contrib_daily_structures_destroyed,
  #               smokePM:ZR, nonsmokePM_MF) %>% 
  # filter_at(vars(AL, AS, BR, CA, CHL, CL, CR, CU,
  #                EC, FE, K, MG, MN, `NA`, NI, NO3, OC, P, PB,
  #                RB, S, SE, SI, SO4, SR, TI, V, ZN, ZR),
  #           any_vars(!is.na(.))) 

future::plan(NULL)

return(merged_smoke_spec_burned_structures_df)
}

# for IVAN
# grid_smoke_structure_dt <- gridcell_smokepm_dt %>%
#   left_join(globfire_structure_joined_df %>%
#               mutate(Id = as.character(Id),
#                      fire_duration=as.integer(FDate-IDate)) %>%
#               mutate(fire_duration = ifelse(fire_duration == 0, 1, fire_duration), # for fires that were only 1 day, calculation gives 0, but should be 1
#                      daily_structures_destroyed=structures_destroyed/fire_duration) %>%
#               dplyr::select(Id, year, IDate, FDate, fire_duration, incid_name_lower, coverage_perc,
#                             structures_destroyed, daily_structures_destroyed),
#             by=c("fire_id"="Id", "year"="year")) %>%
#   mutate(structures_destroyed=ifelse(is.na(structures_destroyed), 0, structures_destroyed),
#          contrib_structures_destroyed=share*structures_destroyed,
#          daily_structures_destroyed=ifelse(is.na(daily_structures_destroyed), 0, daily_structures_destroyed),
#          contrib_daily_structures_destroyed=share*daily_structures_destroyed) %>% 
#   dplyr::select(-smokePM_pred)
# 
# write_feather(grid_smoke_structure_dt, "/Users/ekrasovich/Desktop/grid_smoke_structure_dt.feather")
