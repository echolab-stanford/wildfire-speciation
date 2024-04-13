# Jeff Wen, Emma Krasovich Southworth
# Last Updated: April 2, 2024 by EKS
# Description: processing burned structures dataset

# inputs:
# globfire_fp <- file.path(data_fp, 'intermediate/globfire/globfire_na_final_area_2006-2020.shp')
# mtbs_fp <- file.path(data_fp, 'intermediate/mtbs/mtbs_perims_DD.shp')
# nifc_fp <- file.path(data_fp, "intermediate/fire_locations/WFIGS_-_Wildland_Fire_Locations_Full_History.csv")
# damaged_struct_fp <- file.path(data_fp, 'clean/HE_Structures_Destroyed_2022.xlsx')

processing_burned_structures_data <- function(globfire_fp, mtbs_fp, nifc_fp, damaged_struct_fp) {
  
  # parellelize
  future::plan(multisession) 
  
# shapes ----
## us shape
conus_bbox <- sf::st_bbox(c(xmin=-127.089844, 
                            ymin=22.066785, 
                            xmax=-66.533203, 
                            ymax=50.120578), 
                          crs=st_crs(4326)) %>% 
  st_as_sfc()

us_shape <- tigris::nation(resolution="20m") %>% 
  st_transform(crs=st_crs("epsg:4326")) %>%
  st_crop(conus_bbox)%>%
  st_transform(crs=st_crs("epsg:5070"))

# GlobFire and MTBS loading ----
## load globfire data
globfire_df <- st_read(globfire_fp) %>% 
  st_transform("epsg:5070") %>% 
  st_filter(us_shape) %>% 
  mutate(year=year(IDate)) %>% 
  st_make_valid()

## load mtbs data
mtbs_df <- read_sf(mtbs_fp) %>% 
  filter(Ig_Date>="2006-04-19", Ig_Date<='2020-12-31') %>% 
  st_transform("epsg:5070") %>% 
  mutate(year=year(Ig_Date)) %>% 
  st_make_valid()

# NIFC and structures dataset ----
## read in nifc data to get lat lon for some fires that have integer value lat lon in the structure data
## 515 duplicated unique ids where some are actual duplicates and others are different based on county name, fips, etc.
nifc_fires <- readr::read_csv(nifc_fp) %>% 
  mutate(incident_number=str_replace(UniqueFireIdentifier,"(\\d{4})-(\\w{2})(\\w*)-(\\w*)","\\2-\\3-\\4"),
         year=as.numeric(str_sub(UniqueFireIdentifier, start=1, end=4))) %>% 
  dplyr::select(UniqueFireIdentifier, year, incident_number, IncidentName,
                POOState, POOCounty, POOFips,
                X, Y, InitialLongitude, InitialLatitude,
                TotalIncidentPersonnel, EstimatedCostToDate) %>% 
  arrange(UniqueFireIdentifier) %>% 
  group_by(UniqueFireIdentifier, year, incident_number, 
           POOState, POOCounty, POOFips) %>% 
  summarize(X = mean(X, na.rm=T),
            Y = mean(Y, na.rm=T),
            InitialLongitude = ifelse(all(is.na(InitialLongitude)), NA, mean(InitialLongitude, na.rm=T)),
            InitialLatitude = ifelse(all(is.na(InitialLatitude)), NA, mean(InitialLatitude, na.rm=T)),
            TotalIncidentPersonnel = ifelse(all(is.na(TotalIncidentPersonnel)), NA, max(TotalIncidentPersonnel, na.rm=T)),
            EstimatedCostToDate = ifelse(all(is.na(EstimatedCostToDate)), NA, max(EstimatedCostToDate, na.rm=T))) %>% 
  as.data.frame()

## remove duplicated fires even after summarizing above (about 35 fires with same unique fire ids)
nifc_fires <- nifc_fires %>% 
  filter(UniqueFireIdentifier%not_in%unique(nifc_fires[(duplicated(nifc_fires$UniqueFireIdentifier)) & 
                                                         (!is.na(nifc_fires$UniqueFireIdentifier)),
                                                       "UniqueFireIdentifier"]))

## filter to relevant years and merge in matched lat long from nifc data
## remove 53 fires that still have integer value lat lon
damaged_structure_df <- readxl::read_xlsx(damaged_struct_fp, sheet='Data') %>% 
  as.data.frame() %>% 
  mutate(longitude=-1*longitude) %>% 
  rename(year=yr) %>% 
  filter(year%in%c(2006:2020)) %>%
  left_join(nifc_fires, by=c("incident_number", "year")) %>% 
  mutate(longitude=ifelse((longitude%%1==0) & (!is.na(X)), X, longitude),
         latitude=ifelse((latitude%%1==0) & (!is.na(Y)), Y, latitude)) %>% 
  filter(longitude%%1!=0 & latitude%%1!=0) %>% 
  st_as_sf(coords=c("longitude","latitude"), crs=4326, agr="constant") %>% 
  st_transform("epsg:5070") %>% 
  st_filter(us_shape) %>% 
  mutate(start_date=as.Date(start_date, format="%m/%d/%Y"),
         incident_name_lower=str_to_lower(incident_name))

## spatial join by area first filter to relevant years then only keep obs if enough area coverage
COVERAGE_THRESHOLD <- 0.75
year_list <- 2006:2020

globfire_overlap_df <- purrr::map_df(year_list, function(main_year){
  
  ## subset data to matching year
  temp_globfire_df <- globfire_df %>% 
    filter(year==main_year) %>% 
    mutate(burn_area=st_area(geometry))
  
  temp_mtbs_df <- mtbs_df %>% 
    filter(year==main_year) %>% 
    mutate(mtbs_burn_area=st_area(geometry))
  
  ## calculate the area of intersection 
  temp_intersection_df <- temp_globfire_df %>% 
    st_intersection(temp_mtbs_df) %>% 
    mutate(intersection_area=st_area(geometry)) %>% 
    dplyr::select(Id, Event_ID, irwinID, Incid_Name, Incid_Type, Map_ID,
                  Map_Prog, Asmnt_Type, BurnBndAc, BurnBndLat, BurnBndLon,
                  Ig_Date, Pre_ID, Post_ID, Perim_ID, mtbs_burn_area, 
                  intersection_area) %>% 
    st_drop_geometry() %>% 
    group_by(Id) %>%
    arrange(desc(intersection_area)) %>% 
    filter(row_number()==1) 
  
  ## join to bring in the area of intersection
  temp_globfire_dt <- merge(temp_globfire_df, temp_intersection_df, by="Id", all.x=T) %>% 
    as.data.table()
  
  ## only keep MTBS data if high area coverage overlap and ignition date within start and end globfire date
  SUBSET_COLS <- c("Event_ID","irwinID","Incid_Name","Incid_Type","Map_ID",
                   "Map_Prog","Asmnt_Type","BurnBndAc","BurnBndLat","BurnBndLon",
                   "Ig_Date","Pre_ID","Post_ID","Perim_ID")
  
  ## for obs that have low coverage perc and with fire date outside of the globfire start and end window with additional week boundary on start and end
  ## set values to NA since its probably not a good match
  temp_globfire_dt[,`:=`(coverage_perc=as.numeric(intersection_area/burn_area),
                         mtbs_coverage_perc=as.numeric(intersection_area/mtbs_burn_area))
  ][(coverage_perc<COVERAGE_THRESHOLD) | !((IDate - 7 < Ig_Date) & (FDate + 7 > Ig_Date)),
    (SUBSET_COLS) := lapply(.SD, function(x) NA),
    .SDcols=SUBSET_COLS
  ]
  
}) %>% 
  bind_rows() %>% 
  as.data.frame() %>% 
  st_as_sf()

globfire_structure_joined_df <- globfire_overlap_df %>% 
  mutate(incid_name_lower = str_to_lower(Incid_Name),
         state=str_sub(Event_ID, start=1, end=2)) %>% 
  filter(incid_name_lower %not_in% c(NA, "unnamed")) %>% 
  left_join(damaged_structure_df %>% 
              st_drop_geometry(), 
            by=c("year", "state"="st", "incid_name_lower"="incident_name_lower")) %>% 
  filter(!is.na(structures_destroyed),
         ((IDate-7) <= start_date) & (start_date<=IDate+7)) %>% 
  dplyr::select(Id, year, IDate, FDate, 
                Event_ID, irwinID, incid_name_lower, Incid_Name, Incid_Type, Ig_Date, start_date, 
                BurnBndAc, BurnBndLat, BurnBndLon,Pre_ID, Post_ID, Perim_ID,
                incident_number, incident_name, burn_area, 
                intersection_area, coverage_perc, structures_destroyed) %>% 
  st_drop_geometry()

future::plan(NULL)

## 5 duplicated fires so selected one with more structures destroyed 
globfire_structure_joined_df$Id[duplicated(globfire_structure_joined_df$Id)]

globfire_structure_joined_df <- globfire_structure_joined_df %>% 
  group_by(Id) %>% 
  slice(which.max(structures_destroyed))

return(globfire_structure_joined_df)
}
