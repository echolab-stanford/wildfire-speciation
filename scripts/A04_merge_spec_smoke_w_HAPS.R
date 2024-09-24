# # # Emma Krasovich Southworth, emmars@stanford.edu
# # # Last Updated: April 2, 2024
# # # Description: run regression on burned structures and smoke speciation + plot response curves at different levels of smoke
# 
# function
merge_spec_smoke_w_HAPS <- function(spec_w_smoke_pm_df, haps_spec_df) {

all_spec_df <- spec_w_smoke_pm_df %>%
  left_join(haps_spec_df %>%
              dplyr::select(-Longitude, -Latitude, -epsg, -units) %>%
              mutate(Date = as.Date(Date, format = '%Y-%d-%m')))

return(all_spec_df)

}

