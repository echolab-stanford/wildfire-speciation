# Description to explore data + compare to Kara's

# load the cleaned data
# loadd(pm_spec4data_exp, cache = drake::drake_cache("scripts/.drake"))

# read in state shapefile for mapping -----
us_states <- st_read(file.path(wip_gdrive_fp, 'raw/cb_2018_us_state_500k/cb_2018_us_state_500k.shp')) %>% 
  st_transform(4326) %>% 
  dplyr::select(State = STUSPS, st_name = NAME)


# how long are monitors monitoring?
num_yrs_site_data <- pm_spec4data_exp %>% 
  group_by(SiteCode, Dataset, State, Latitude, Longitude) %>% 
  dplyr::summarise(min_yr = min(year),
                   max_yr = max(year),
                   duration = (max_yr - min_yr) + 1) %>% 
  ungroup() %>% 
  mutate(yr_cat = case_when(
    duration >= 15 ~ '15 years, full sample',
    duration >= 10 & duration <= 14 ~ '10-14 years',
    duration >= 5 & duration <= 9 ~ '5-9 years',
    duration < 5 ~ '< 5',
  )) %>% 
  st_as_sf(coords =c(x = 'Longitude', y = 'Latitude'), crs = 4326)


num_yrs_site_data$yr_cat <- factor(num_yrs_site_data$yr_cat, 
                                   order = TRUE, 
                                   levels =c('< 5', '5-9 years', '10-14 years', '15 years, full sample'))


mapview(num_yrs_site_data, zcol= 'yr_cat', cex = 5, 
        col.regions=list("yellow","orange", "darkred", "darkslateblue"))



#----------------------------------------------------------------------------------------------
# One test that might be worth trying would be to calculate the % of components/MF
# (e.g., SO4/MF and OC/MF), and see if they add up to close to one




# find what fraction of species are negative ------------------------------------------------
neg_spec_pct <- pm_spec4data_exp %>% 
  dplyr::select(Dataset, Date, State, region, season, year, smokePM_pred, MF, RC_PM2, 
                NH4, tot_metals, ammNO3, ammSO4, AL, AS, BR, CA, EC, OC, 
                CHL, CL, CR, CU, FE, K, MG, MN, MO, `NA`, NI, NO3, N2, P, 
                PB, RB, S, SE, SI, SOIL, SO4, SR, V, ZN, ZR) 
  
var_list <- c('smokePM_pred', 'MF', 'RC_PM2', 'NH4', 'tot_metals', 'ammNO3', 'ammSO4', 'AL', 
              'AS', 'BR', 'CA', 'EC', 'OC', 'CHL', 'CL', 'CR', 'CU', 
              'FE', 'K', 'MG', 'MN', 'MO', 'NA', 'NI', 'NO3', 'N2', 'P', 
              'PB', 'RB', 'S', 'SE', 'SI', 'SOIL', 'SO4', 'SR', 'V', 'ZN', 'ZR')

# current_var <- var_list[1] # test one stat

summ_stats <- purrr::map_df(var_list, function(current_var) {
  
 current_df <- neg_spec_pct[names(neg_spec_pct) == current_var] 

  current_stats <- tibble(
    var = current_var,
    num_missing = sum(is.na(current_df)),
    num_neg = length(which(current_df < 0)),
    tot_obs_not_missing = nrow(current_df) - sum(is.na(current_df))
  ) 
  
}) %>% 
  bind_rows() %>% 
  mutate(pct_neg_of_nonNA = round((num_neg/tot_obs_not_missing)*100, digits = 2),
         pct_missing = round((num_missing/(tot_obs_not_missing + num_missing))*100, digits = 2)) %>% 
  mutate(neg_cat = as.factor(ifelse(pct_neg_of_nonNA >= 20, '>= 20% negative', '< 20% negative'))) %>% 
  mutate(missing_cat = as.factor(ifelse(pct_missing >= 20, '>= 20% missing', '< 20% missing'))) %>% 
  arrange(var)
  

# make plot of how many negatives in data that is non missing
pct_spec_neg_plot <- ggplot(summ_stats,
                       aes(x = var, y = pct_neg_of_nonNA, fill = 'steelblue')) +
  scale_fill_manual(values=c("steelblue")) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 10, color = 'red') +
  geom_text(aes(label= round(pct_neg_of_nonNA, digits = 1)), vjust=-1) +
  xlab('Species') + 
  ylab('Percent (%)') +
  ggtitle('% of reported, non-missing species data that is negative') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
pct_spec_neg_plot

# make plot of how much data
pct_spec_missing_plot <- ggplot(summ_stats,
                            aes(x = var, y = pct_missing, fill = 'steelblue')) +
  scale_fill_manual(values=c("forestgreen",
                             "steelblue")) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 50, color = 'red') +
  geom_text(aes(label= round(pct_missing, digits = 1)), vjust=-1) +
  xlab('Species') + 
  ylab('Percent (%)') +
  ggtitle('% of species data that is missing for all dates in sample') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
pct_spec_missing_plot

# -------------------------------------------------------------------------------------
# plot negatives by site

# how long are monitors monitoring?
full_sample_site_avgs <- pm_spec4data_exp %>% 
  group_by(SiteCode, Dataset, State, Latitude, Longitude) %>% 
  dplyr::summarise(min_yr = min(year),
                   max_yr = max(year),
                   duration = (max_yr - min_yr) + 1) %>% 
  ungroup() %>% 
  mutate(yr_cat = case_when(
    duration >= 15 ~ '15 years, full sample',
    duration >= 10 & duration <= 14 ~ '10-14 years',
    duration >= 5 & duration <= 9 ~ '5-9 years',
    duration < 5 ~ '< 5',
  )) %>% 
  st_as_sf(coords =c(x = 'Longitude', y = 'Latitude'), crs = 4326)









# find the difference between reconsituted and measured PM2.5 mass -------------
pm25_site_averages_year <- pm_spec4data_exp %>% 
  group_by(Dataset, SiteCode, Latitude, Longitude, epsg, year) %>% 
  dplyr::summarise(MF_yr_avg = mean(MF, na.rm = TRUE),
                   RC_PM2_yr_avg = mean(RC_PM2, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(PM25_ratio = RC_PM2_yr_avg/MF_yr_avg) %>% 
  # mutate(date = as.Date(paste(year, month, 1, sep = "-"))) %>% 
  filter(!is.nan(PM25_ratio))

# randomly sample 10 sites for display purposes
set.seed(1)
site_sample <- sample(pm25_site_averages$SiteCode, 20)

pm25_diff_ts <- ggplot(pm25_site_averages %>% 
                         filter(SiteCode %in% site_sample), 
                       aes(x = year, y = PM25_ratio)) +
  geom_point() +
  xlab('year') + ylab('Ratio bw Recon + Measured PM2.5 for sample of sites') +
  facet_wrap(~SiteCode) +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  theme_bw()
pm25_diff_ts

# what about when we use Marissa's measure of smoke pm
cont_pm_spec_fracs <- pm_spec4data_exp %>% 
  mutate(OC_frac = ifelse(smokePM_pred != 0, OC_adj/smokePM_pred, NA),
         EC_frac = ifelse(smokePM_pred != 0, EC_adj/smokePM_pred, NA), 
         NO3_frac = ifelse(smokePM_pred != 0, NO3_adj/smokePM_pred, NA), 
         dust_frac = ifelse(smokePM_pred != 0, SOIL_adj/smokePM_pred, NA), 
         SO4_frac = ifelse(smokePM_pred != 0, SO4_adj/smokePM_pred, NA), 
         NH4_frac = ifelse(smokePM_pred != 0, NH4/smokePM_pred, NA), 
         metals_frac = ifelse(smokePM_pred != 0, tot_metals/smokePM_pred, NA)) 

test <- cont_pm_spec_fracs %>% 
  dplyr::select(Dataset, SiteCode, State, region, season, year, month,
                smokePM_pred, RC_PM2, MF_adj, OC_frac, EC_frac, NO3_frac, 
                dust_frac, SO4_frac, NH4_frac, metals_frac)
               
                             

#Descriptive stats:
# get monthly averages for categories we care about:
monthly_site_species <- pm_spec4data_exp %>% 
  dplyr::select(Dataset, SiteCode, State, region, season, year, month, smoke_day, 
                low_count, med_count, high_count, density_missing, smokePM_pred, RC_PM2,
                NH4, tot_metals, OC, EC, NO3, SOIL, SO4, ammNO3, ammSO4, Latitude, Longitude, epsg) %>% 
  group_by(Dataset, SiteCode, year) %>% 
  mutate(num_site_smoke_days_yrly = sum(smoke_day == 1, na.rm = TRUE),
         num_site_nonsmoke_days_yrly = sum(smoke_day == 0, na.rm = TRUE)) %>% 
  mutate(smoke2nonsmokeday_yrly_ratio = num_site_smoke_days_yrly/num_site_nonsmoke_days_yrly) %>% 
  ungroup() %>% 
  # create monthly avg of species
  group_by(Dataset, SiteCode, season, year, month, ) %>% 
  group_by(region, State, Dataset, SiteCode, season, year, month, Latitude, Longitude, epsg)
