# Emma Krasovich Southworth, Marissa Childs
# 11/19/24
# read in Childs et al (2022) data on CSN and IMPROVE
# Creates a data frame of daily EPA station total PM2.5 (2000-2022) and smoke 

# ------------------------------------------------------------------------------
# # R1, Q1: COMPARING CHILDS ET AL 2022 TO STATION BASED MEASURE OF SMOKE
# -----------------------------------------------------------------------------
# REFERENCE FUNCTIONS -----------------------------------------------------
# set up function from childs et al (2022)
nonsmoke_medians <- function(data, 
                             main_var, 
                             smoke_day, 
                             spatial_unit, 
                             temporal_unit, 
                             temporal_trend){
  main_var <- enquo(main_var)
  smoke_day <- enquo(smoke_day)
  spatial_unit <- enquo(spatial_unit)
  temporal_unit <- enquo(temporal_unit)
  temporal_trend <- enquo(temporal_trend)
  
  new_name <- paste0(rlang::as_name(main_var), "_med_3yr")
  
  full_panel <- expand.grid(site_id = data %>% pull(!!spatial_unit) %>% unique, 
                            month = data %>% pull(!!temporal_unit) %>% unique, 
                            year = data %>% pull(!!temporal_trend) %>% unique) %>% 
    rename(!!spatial_unit := site_id, 
           !!temporal_unit := month, 
           !!temporal_trend := year) %>% 
    ungroup
  
  data %>% 
    filter(!is.na(!!main_var) & !!smoke_day == 0) %>%
    full_join(full_panel) %>%
    group_by(!!spatial_unit, !!temporal_unit, !!temporal_trend) %>%
    summarise(main_var = list(!!main_var),
              nobs = n(),
              .groups = "drop") %>%
    arrange(!!spatial_unit, !!temporal_unit, !!temporal_trend) %>%
    group_by(!!spatial_unit, !!temporal_unit) %>%
    mutate(main_var_lag = lag(main_var, n = 1, default = list(NA)),
           main_var_lead = lead(main_var, n = 1, default = list(NA)),
           nobs_lag = lag(nobs, n = 1, default = 0),
           nobs_lead = lead(nobs, n = 1, default = 0)) %>%
    ungroup %>%
    rowwise %>%
    mutate(main_var_3yr = list(c(main_var, main_var_lag, main_var_lead)),
           main_var_med_3yr = median(unlist(main_var_3yr), na.rm = T),
           nobs_3yr = nobs + nobs_lead + nobs_lag) %>%
    ungroup %>%
    transmute(!!spatial_unit, !!temporal_unit, !!temporal_trend,
              nobs_3yr,
              !!new_name := main_var_med_3yr)
}

# CALCULATE THE STATION SMOKE DATA FROM OUR DATASET  ----------------------------
# Load PM2.5 from EPA stations
epa_pm <- read.csv(file.path(data_fp, 'clean/station_PMspeciation_2006_2020.csv')) %>% 
  dplyr::select(Dataset:month, monitor_month:epsg) %>% 
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>% 
  dplyr::select(state_name:region, year, month, monitor_month:site_id, pm25= 'MF_adj')

# # Load station locations
epa_ll <- read.csv(file.path(data_fp, 'clean/final_PMspeciation_sites.csv')) %>% 
  st_as_sf(coords = c('long', 'lat'), crs = 4326) %>% 
  st_make_valid()

# Set sample range
sample_start = ymd("2006-01-01")
sample_end = ymd("2021-01-01")

# Expand to full station-day panel
epa_panel = expand.grid(site_id = epa_ll$site_id, 
                        Date = seq.Date(sample_start, sample_end, by = "day")) %>% 
  left_join(epa_pm, by = c("site_id", "Date"))

# Load smoke data
smoke = readRDS(file.path(data_fp, "revisions/smoke_plumes_sfdf_20060101-20230705.RDS")) %>% 
  st_transform(4326) %>% 
  st_make_valid() %>% 
  filter(st_is_valid(.)) # & between(st_bbox(.)$xmin, -180, 180) & between(st_bbox(.)$ymin, -90, 90)) %>% 
 

smoke_missing_dates = readRDS(file.path(data_fp, "revisions/smoke_dates_not_online_20060101-20230705.RDS")) 
plume_start = min(ymd(unique(smoke$date)))
plume_end = max(ymd(unique(smoke$date)))

# Get smoke days at station locations
# Takes 1-2 minutes
station_smoke <- st_intersects(epa_ll, smoke) %>% 
  map(function(i) {
    smoke[i,] %>% pull(date) %>% unique()
  })

station_smokedays <- epa_ll %>% 
  st_drop_geometry() %>%
  select(site_id) %>% 
  mutate(Date = station_smoke) %>% 
  unnest(Date) %>% 
  mutate(Date = ymd(Date), 
         plume = 1) %>% 
  as.data.frame()

# Calculate smoke PM2.5
smoke_pm <- epa_panel %>% 
  mutate(year = year(Date), month = month(Date)) %>% 
  left_join(station_smokedays, by = c('site_id', 'Date')) %>% 
  # Replace missing plume info with 0s, except on missing smoke data days, or 
  # outside of plume date ranges
  replace_na(list(plume = 0)) %>% 
  mutate(plume = ifelse(Date %in% smoke_missing_dates, NA, plume), 
         plume = ifelse((Date > plume_end) | (Date < plume_start), NA, plume)) %>% 
  # Calculate non-smoke medians, i.e. background PM
  # Drop station-days not observed, dates before start of plume data, and dates 
  # where plume data were otherwise not online
  {left_join(., 
             nonsmoke_medians(filter(., !is.na(pm25), !is.na(plume)), 
                              pm25, plume, site_id, month, year), 
             by = c("site_id", "month", "year"))} %>% 
  # Fill in smoke days
  mutate(smoke_day = ifelse(Date %in% smoke_missing_dates, 0, plume), 
         plume_lag1 = lag(plume, 1), 
         plume_lag2 = lag(plume, 2), 
         plume_lead1 = lead(plume, 1), 
         plume_lead2 = lead(plume, 2), 
         plume_tnn = ifelse(Date %in% smoke_missing_dates, 
                            as.integer(ifelse(
                              (!is.na(plume_lag1)) & (!is.na(plume_lead1)),
                              (plume_lag1 == 1) & (plume_lead1 == 1),
                              ifelse(is.na(plume_lag1), T, plume_lag1 == 1) & 
                                ifelse(is.na(plume_lead1), T, plume_lead1 == 1) & 
                                ifelse(is.na(plume_lag2), T, plume_lag2 == 1) & 
                                ifelse(is.na(plume_lead2), T, plume_lead2 == 1)
                            )), 
                            plume), 
         filled_smoke_day = plume_tnn, 
         # Calculate PM anomalies and smokePM
         pm25_anom = pm25 - pm25_med_3yr, 
         smokePM = ifelse(smoke_day, pmax(pm25_anom, 0), 0), 
         filled_smokePM = ifelse(filled_smoke_day, pmax(pm25_anom, 0), 0)) %>% 
  left_join(epa_ll %>% 
              st_drop_geometry(), 
            by = "site_id") 

# Save
# saveRDS(smoke_pm, 
#         file.path(data_fp, "revisions/estimated_station_smokePM_eks.rds")
# )
# drop to save memory
rm(station_smoke, smoke, smoke_pm)

# COMBINE OUR STATION BASED MEASURE OF SMOKE WITH PREDICTED VALUES ----------------------------
# station PM
station_df <- readRDS(file.path(data_fp, '/revisions/estimated_station_smokePM_eks.rds')) %>% 
  dplyr::select(site_id, Date, year, month, monitor_month, pm25, smokePM)
# need to merge with the grid cell id

# read in predictions
preds_df <- readRDS(file.path(data_fp, 'revisions/smokePM_predictions_20060101_20201231.rds')) %>% 
  rename(Date = 'date')

# Load 10 km grid
grid_10km <- st_read(file.path(data_fp, 'intermediate/grid/grid_10km_wgs84.shp')) %>% 
  st_transform(4326)

# Match monitor locations in speciation data with grid 10km ID ####
# pull id cells that intersect with the monitoring sites
pts_in_grid_df <- epa_ll %>% 
  st_join(grid_10km) %>% 
  st_drop_geometry() %>% 
  rename(grid_id_10km = 'ID') %>% 
  dplyr::select(-epsg, -COORDX, -COORDY)

# Add grid_id_10km to the speciation df, join the smoke PM data by grid cell
station_gridded_df <- station_df %>% 
  left_join(pts_in_grid_df, 
            by = c("site_id")) %>% 
  # match PM2.5 from the smoke_pm dataset by gridcell
  left_join(preds_df, 
            by = c("grid_id_10km", "Date")) %>% 
  # Replace NA values (non-smoke days) to 0 (which would mean no smoke detected)
  mutate(smokePM_pred = ifelse(is.na(smokePM_pred), 0, smokePM_pred)) %>% 
  distinct() %>% 
  dplyr::select(site_id:pm25, station_smokePM = 'smokePM', grid_id_10km, smokePM_pred)


# RUN REGRESSION AND CORRELATION TO CHECK RELATIONSHIP -------------------------
comp_smoke_reg <- feols(smokePM_pred ~ station_smokePM | site_id,
                        data = station_gridded_df)
# etable(comp_smoke_reg)
summary(comp_smoke_reg)

# PLOT IN SCATTER!!!!  ---------------------------------------------------------
# Create the scatter plot
scatter_plot <- ggplot(station_gridded_df, 
                       aes(x = station_smokePM, y = smokePM_pred)) +
  geom_point(size = 2, alpha = 0.7) +  # Scatter plot points with size and transparency
  labs(
    x = "Constructed smoke PM2.5 for stations in our sample",
    y = "Predicted Smoke PM2.5 in Childs et al (2022) ") +
  theme_minimal() +  # Minimalistic theme
  geom_abline(slope = 1, intercept = 0, color = "coral", linetype = "dashed") +  # Add y = x line
  theme(plot.title = element_text(hjust = 0.5), # Center align the title
    legend.position = "right"              # Position of the legend
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18), # Center align and enlarge title
    axis.title = element_text(size = 16),             # Increase size of axis titles
    axis.text = element_text(size = 14),              # Increase size of axis text
    legend.text = element_text(size = 12),            # Increase size of legend text
    legend.title = element_text(size = 14),           # Increase size of legend title
    legend.position = "right"                         # Position of the legend
  ) +
  ggpubr::stat_cor(
    aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
    method = "pearson",
    label.x = 5, label.y = 400, # Position of correlation label
    size = 6
  )  # Add correlation
scatter_plot

# save file
ggsave(
  filename = 'SIFig3_station_vs_preds_smoke_corr.png',
  plot = scatter_plot,
  path = file.path(results_fp, 'SI Figs'),
  scale = 1,
  width = 8,
  height = 8,
  dpi = 320)

# save file
# ggsave(
#   filename = 'SIFig3_station_preds_smoke_corr.pdf',
#   plot = scatter_plot,
#   path = file.path(results_fp, 'SI Figs'),
#   scale = 1,
#   width = 8,
#   height = 8)

# COMBINE OUR SUBSET OF STATIONS WITH MARISSAS SUBSET OF STATIONS ---------------

# stations
# childs <- readRDS(file.path(data_fp, "revisions/station_smokePM.rds")) 
# 
# combined_check_stations <- station_gridded_df %>% 
#   left_join(childs) %>% 
#   distinct(site_id, Date,grid_id_10km, station_smokePM, smokePM)
# 
# # Create the scatter plot
# scatter_plot2 <- ggplot(combined_check_stations, 
#                        aes(x = station_smokePM, y = smokePM)) +
#   geom_point(size = 2, alpha = 0.7) +  # Scatter plot points with size and transparency
#   labs(
#     x = "Constructed smoke PM2.5 for stations in our sample",
#     y = "Constructed station smoke PM2.5 in Childs et al (2022)") +
#   theme_minimal() +  # Minimalistic theme
#   geom_abline(slope = 1, intercept = 0, color = "coral", linetype = "dashed") +  # Add y = x line
#   theme(plot.title = element_text(hjust = 0.5), # Center align the title
#         legend.position = "right"              # Position of the legend
#   ) #+
#   # ggpubr::stat_cor(
#   #   aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
#   #   method = "pearson",
#   #   label.x = 0, label.y = 100, # Position of correlation label
#   #   size = 5
#   # )  # Add correlation
# scatter_plot2
# 
# # save file
# ggsave(
#   filename = 'SIFig_station_smoke_comparison_corr_eks.png',
#   plot = scatter_plot2,
#   path = file.path(results_fp, 'SI Figs'),
#   scale = 1,
#   width = 12,
#   height = 8,
#   dpi = 320)


# ----------------------------------------------------------------------------
# R1, Q2: VALIDATE OUR MODEL BY INCLUDING ADDITIONAL COVARIATES
# ----------------------------------------------------------------------------
# read in precipation data from CHIRPS
ppt_files <- list.files(file.path(gee_data_fp), pattern = 'precip', full.names = TRUE)
# current_fp <- ppt_files[1]

ppt_df <- map_df(ppt_files, function(current_fp) {
  
  current_df <- read.csv(current_fp) %>% 
    dplyr::select(grid_id_10km = 'ID', ppt_mm_per_day = 'sum') %>% 
    mutate(year = strapply(current_fp, "\\d{4}", as.numeric, simplify = TRUE))
  
})

# read in temperature data from ERA5
temp_files <- list.files(file.path(gee_data_fp), pattern = 'temp', full.names = TRUE)
# current_fp <- temp_files[1]

temp_df <- map_df(temp_files, function(current_fp) {
  
  current_df <- read.csv(current_fp) %>% 
    dplyr::select(grid_id_10km = 'ID', temp = 'mean_2m_air_temperature') %>% 
    mutate(year = strapply(current_fp, "\\d{4}", as.numeric, simplify = TRUE)) %>% 
    mutate(tempF = kelvin.to.fahrenheit(temp, round = 2))

})
 
covar_df <-  left_join(temp_df, ppt_df) %>% 
  dplyr::select(-temp)
  
# load these datasets
loadd(c(clean_PMspec_df,clean_PMspec_sites_df, parameter_categories), cache = drake_cache)

# Load 10 km grid
grid_10km <- st_read(file.path(data_fp, 'intermediate/grid/grid_10km_wgs84.shp')) %>% 
  st_transform(4326)

# Match monitor locations in speciation data with grid 10km ID ####
# pull id cells that intersect with the monitoring sites
pts_in_grid_df <- clean_PMspec_sites_df %>% 
  st_as_sf(coords = c('long', 'lat'), crs = 4326) %>% 
  st_join(grid_10km) %>% 
  st_drop_geometry() %>% 
  rename(grid_id_10km = 'ID') %>% 
  dplyr::select(-epsg, -COORDX, -COORDY)

# Add grid_id_10km to the speciation df, join the smoke PM data by grid cell
station_gridded_df <- clean_PMspec_df %>% 
  left_join(pts_in_grid_df, 
            by = c("site_id", 'Dataset')) %>% 
  left_join(covar_df) %>% 
  dplyr::select(Dataset:region, year, month, Date, 
                monitor_month, grid_id_10km, smoke_day, site_id, 
                tempF, ppt_mm_per_day, MF_adj, smokePM, 
                nonsmokePM_MF, AL:ZN) %>% 
  # drop Soil and zirconium, not needed for this analysis
  dplyr::select(-SOIL)

# GET BASELINE AVERAGES FOR EACH SPECIES --------------------------------
# calculate the baseline for each chemical species across full sample by filtering to nonsmoke days
# rather than all days, bc some places may be really smoky 
# relative to others 
spec_ns_samp_avgs_df <- station_gridded_df %>% 
  dplyr::select(site_id, Date, AL:ZN, smoke_day) %>% 
  pivot_longer(cols = c(AL:ZN), names_to ='species', values_to = 'conc_val') %>% 
  filter(smoke_day == 'nonsmoke day') %>% 
  group_by(species) %>% 
  dplyr::summarise(avg_nonsmoke_spec_conc = mean(conc_val, na.rm = TRUE), .groups = 'drop')


# RUN REGRESSION FOR SPECIES -------------------------------------------------
# RUN IN LEVELS ACROSS FULL SAMPLE, 
# DIVIDE BETAS BY SAMPLE AVG TO GET % CHANGE RELATIVE TO SAMPLE BASELINE

covar_reg = feols(c(AL,AS,BR, CA, CL, CR, CU, EC, FE, 
                            K, MG,MN, `NA`, NI, NO3, OC, P,  PB, RB,
                            S,  SE, SI, SO4, SR, TI, V,  ZN)
                          ~ smokePM + nonsmokePM_MF + tempF + ppt_mm_per_day | 
                            monitor_month + year, station_gridded_df, cluster = 'site_id')
# summary(covar_reg) # look at regression results
## calculate 95 CI%
CIs <- confint(
  covar_reg) %>% 
  rename(species = 'lhs',
         pm_type = 'coefficient',
         CI25 = '2.5 %',
         CI975 = '97.5 %') %>% 
  mutate(measure = 'MF') %>% 
  dplyr::select(-id)

# get coefficients and prepare for plotting
covar_coeffs <- coeftable(covar_reg) %>% 
  rename(pval = 'Pr(>|t|)',
         se = 'Std. Error',
         species = 'lhs',
         pm_type = 'coefficient') %>% 
  mutate(pval = round(pval, digits = 4)) %>% 
  # get pvalues
  mutate(sig = ifelse(pval > .05, 'non-significant', 'significant')) %>% 
  dplyr::select(-id)  %>% 
  left_join(parameter_categories, by = 'species') %>% 
  mutate(species_type = fct_relevel(species_type,
                                    c("Alkaline-earth metals", "Alkali metals", 
                                      "Transition metals", "Metalloids", "Other metals", 
                                      "Nonmetals",  "Halogens", "Organics"))
  ) %>% 
  mutate(measure = 'MF') %>% 
  left_join(CIs, by = c('species', 'pm_type', 'measure'))


# merge sample avg for each species and divide each species' betas by full sample avg for each species
covar_coeffs_normalized <- covar_coeffs %>% 
  left_join(spec_ns_samp_avgs_df, by = 'species') %>% 
  mutate(norm_est = Estimate/avg_nonsmoke_spec_conc, # how much a species has changed relative to its baseline
         norm_CI25 = CI25/avg_nonsmoke_spec_conc,
         norm_CI975 = CI975/avg_nonsmoke_spec_conc) %>% 
  filter(pm_type == 'smokePM') %>% 
  mutate(model_type = 'environmental covariates included')

# load the original coefficients
loadd(spec_pal, full_samp_PMcoeffs_normalized)

og_coeffs <- full_samp_PMcoeffs_normalized %>% 
  mutate(model_type = 'no environmental covariates')

# add together:
all_coeffs <- bind_rows(covar_coeffs_normalized, og_coeffs)

# plot percent change for all species ------------------------------------------
pct_change_samp_reg_plot <- ggplot(all_coeffs, 
                                   aes(x = species_long,
                                       y = 100*norm_est, 
                                       color=species_type,
                                       shape = model_type
                                   )) +
  geom_point(size=3, alpha = 0.6, stat = "identity") +
  geom_linerange(aes(ymin = (100*norm_CI25), 
                     ymax = (100*norm_CI975)), stat = "identity") +
  scale_x_discrete(limits = c(
    # "Organics"
    "Organic Carbon (OC)", "Elemental Carbon (EC)",
    # "Halogens"
    "Bromine (Br)", "Chlorine (Cl)", 
    #  "Nonmetals"
    "Phosphorus (P)","Sulfur (S)", "Sulfate (SO4)",  "Nitrate (NO3)", "Selenium (Se)", 
    # "Other metals"
    "Aluminum (Al)", "Lead (Pb)",
    # "Metalloids"
    "Silicon (Si)", "Arsenic (As)",
    # "Transition metals"
    "Manganese (Mn)", "Zinc (Zn)", "Titanium (Ti)","Iron (Fe)", "Copper (Cu)", "Vanadium (V)", "Nickel (Ni)", "Chromium (Cr)",
    # "Alkali metals"
    "Potassium (K)", "Rubidium (Rb)", "Sodium (Na)",
    # "Alkaline-earth metals"
    "Strontium (Sr)","Calcium (Ca)",  "Magnesium (Mg)" 
  )) +
  scale_color_manual(values= spec_pal) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey') +
  labs(y = expression(paste('% change relative to average nonsmoke day concentration')),
       x = 'Species',
       shape = 'Model inclusions',
       title = expression(paste("% change in concentration for 1 ug/", m^3, " increase in smoke ", "PM"["2.5"]))) + 
  theme_light() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        title= element_text(size=12, face='bold'),
        axis.title.x = element_text(size=11, face = 'plain'),
        axis.title.y = element_text(size=11, face = 'plain')) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'grey50') +
  guides(color = 'none') 
pct_change_samp_reg_plot


# save file
ggsave(
  filename = 'SIFIG_environmental_covars_vs_og_coeffs.pdf',
  plot = pct_change_samp_reg_plot,
  path = file.path(results_fp, 'SI Figs'),
  scale = 1,
  width = 7,
  height = 10,
  dpi = 320) 

ggsave(
  filename = 'SIFIG_environmental_covars_vs_og_coeffs.png',
  plot = pct_change_samp_reg_plot,
  path = file.path(results_fp, 'SI Figs'),
  scale = 1,
  width = 7,
  height = 10,
  dpi = 320) 



# ---------------------------------------------------------------------
# PLOT ONE CHEMICAL PREDICTE CONC AS AN EXAMPLE:
# ---------------------------------------------------------------------
loadd(attributable_frac_df)
datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-01-01"), by = "24 month")


as_df <- attributable_frac_df %>% 
  filter(species == 'AS') %>% 
  filter(label == 'Attributable concentration to smoke PM2.5')

avg_AS_pred_spec_plot <- ggplot(as_df) +
  geom_line(aes(x = mon_yr,
                   y = conc), 
                   color = "#662D91", size =.6, linetype = 'longdash') +
  labs(title = expression("Arsenic concentration attributable to wildfire smoke PM2.5"),
       y = expression("Concentration (\u00b5g/m"^3*")"),
       x = "") +
  scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
  geom_hline(aes(yintercept = 0)) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0),
        axis.text = element_text(),
        axis.title = element_text(),
        legend.text = element_text(),
        legend.title = element_blank(),
        axis.ticks.x = element_line(linetype = 1),
        axis.ticks.length=unit(-0.2, "cm"),
        panel.grid = element_blank(),
        legend.position = "bottom",
        panel.grid.major.y = element_line(color = "gray95"),
        axis.line.y = element_line(color = "grey10"))
avg_AS_pred_spec_plot

# save file
ggsave(
  filename = 'SIFig1_OUTPUTS_smoke_attributable_predictions_AS.pdf',
  plot = avg_AS_pred_spec_plot,
  path = file.path(results_fp, 'SI Figs'),
  scale = 1,
  width = 6,
  height = 4,
  dpi = 320)



