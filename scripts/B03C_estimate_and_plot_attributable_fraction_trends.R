# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 6, 2023
# Description: make a time series plot for each species

# loadd(c(pred_obs_joined_df, pm_pal, spec_pal, pred_conc_fp, parameter_categories, clean_PMspec_df), cache = drake::drake_cache(".drake"))
# pred_conc_fp <- file.path(data_fp, "intermediate/predicted_concentrations_all_days.fst")

estimate_and_plot_attributable_fraction_trends <- function(pred_obs_joined_df, 
                                                           pm_pal, spec_pal,
                                                           pred_conc_fp, 
                                                           parameter_categories,
                                                           clean_PMspec_df) {

  # ATTRIBUTABLE FRACTION ----------------------------------------------------------
  # calc fraction and determine if the fraction of species concentration attributable to wildfire
  # is significantly increasing over time
  calc_frac <- pred_obs_joined_df %>% 
    mutate(att_frac = pred_spec_smoke/obs_conc) %>% 
    dplyr::select(Date, species, region, year, month, att_frac, site_id, grid_id_10km, monitor_month) %>% 
    pivot_wider(names_from = 'species', values_from = 'att_frac') 
  
  # TEST IF THE ATTRIBUTABLE FRACTION IS INCREASING OVER TIME
  print(paste('running model'))
  att_frac_model = feols(c(AL,AS,BR,CA,CL, CR,CU, EC, FE, K, MG, MN, 
                           `NA`, NI, NO3, OC, P, PB, RB, S, 
                           SE, SI, SO4, SR, TI, V, ZN)
                         ~  year | monitor_month,
                         data = calc_frac, cluster = 'site_id')
  # etable(att_frac_model)
  
  ## calculate 95 CI%
  CIs <- confint(att_frac_model) %>% 
    rename(CI25 = `2.5 %`,
           CI975 = `97.5 %`,
           species = 'lhs') %>% 
    dplyr::select(-id)
  
  # get coefficients and prepare for plotting
  coeffs <- coeftable(att_frac_model) %>%
    rename(pval = 'Pr(>|t|)',
           se = 'Std. Error',
           species = 'lhs',
           coeff = 'Estimate') %>%
    mutate(pval = round(pval, digits = 3)) %>%
    dplyr::select(-id) %>%
    left_join(CIs, by = c('species', 'coefficient')) %>% 
    mutate(sig = ifelse(pval < .05, "sig increasing", "cannot detect increase"))
  
  rm(att_frac_model)
  
  print(paste('averaging preds for each region'))
  
  # combine all days with predicted dataset
  avg_conc_site_month_yr <- read.fst(pred_conc_fp) %>%
    # calculate avg attributable fraction of species concentration to total concentration
    # get avg concentration across all sites for a given month and year
    group_by(region, species, year, month, site_id) %>%
    dplyr::summarise(avg_smoke_conc = mean(pred_spec_smoke, na.rm = TRUE),
                     avg_nonsmoke_conc = mean(pred_spec_nonsmoke, na.rm = TRUE),
                     avg_tot_conc = mean(pred_spec_allPM, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(species, year, month) %>%
    dplyr::summarise(avg_smoke_conc = mean(avg_smoke_conc, na.rm = TRUE),
                     avg_nonsmoke_conc = mean(avg_nonsmoke_conc, na.rm = TRUE),
                     avg_tot_conc = mean(avg_tot_conc, na.rm = TRUE)) %>%
    ungroup()
  
 
  # get monthly avg conc for the observed "total concentration" of a species
  obs_improve_monthly_mean <- clean_PMspec_df %>%
    # select the vars that are needed for the
    dplyr::select(Dataset, region, year, month, Date,
                  smoke_day, site_id, MF_adj, smokePM,
                  nonsmokePM_MF, AL:ZN) %>%
    mutate(monitor_month = paste0(site_id,"_",month)) %>%
    pivot_longer(cols = c(AL:ZN), names_to ='species', values_to = 'conc_val') %>%
    group_by(region, species, year, month, site_id) %>%
    dplyr::summarise(avg_obs_tot_conc = mean(conc_val, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(region,species, year, month) %>%
    dplyr::summarise(avg_obs_tot_conc = mean(avg_obs_tot_conc, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(species, year, month) %>%
    dplyr::summarise(avg_obs_tot_conc = mean(avg_obs_tot_conc, na.rm = TRUE)) %>%
    ungroup()
  
  # now merge back with predicted data
  print(paste('merging predictions with observed data'))
  avg_attributable_fraction_vs_obs <- avg_conc_site_month_yr %>%
    left_join(obs_improve_monthly_mean, 
              by = c('species', 'year', 'month'))
  
  # # prep df for plotting
  plot_df <- avg_attributable_fraction_vs_obs %>%
    dplyr::select(-avg_obs_tot_conc, -avg_nonsmoke_conc, -avg_tot_conc) %>%
    mutate(label = 'Attributable concentration to smoke PM2.5') %>%
    rename(conc = 'avg_smoke_conc') %>%
    bind_rows(avg_attributable_fraction_vs_obs %>%
                dplyr::select(conc = 'avg_obs_tot_conc', year, month, species) %>%
                mutate(label = 'Observed Concentration in Total PM2.5')) %>%
    mutate(mon_yr = as.Date(paste0(year, "-", month, "-", "01"), format = "%Y-%m-%d")) %>%
    left_join(parameter_categories, by = 'species') %>%
    mutate(label = fct_relevel(label, c("Observed Concentration in Total PM2.5",
                                        "Attributable concentration to smoke PM2.5"))) %>%
    filter(!is.na(species)) %>% 
    left_join(coeffs %>% rename(year_trend = 'coeff'), by = 'species')
  
  datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-12-01"), by = "24 month")
  
  limits = c(
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
    "Manganese (Mn)", "Zinc (Zn)", "Titanium (Ti)", "Iron (Fe)", "Copper (Cu)", "Vanadium (V)", "Nickel (Ni)", "Chromium (Cr)",
    # "Alkali metals"
    "Potassium (K)", "Rubidium (Rb)", "Sodium (Na)",
    # "Alkaline-earth metals"
    "Strontium (Sr)","Calcium (Ca)",  "Magnesium (Mg)")
  
  # PLOTTING
  avg_mon_pred_spec_plot <- ggplot(plot_df) +
    geom_density(aes(x = mon_yr,
                     y = conc,
                     fill = label,
                     color = label),
                 stat = "identity", alpha = 0.8) +
    labs(title = expression("Chemical species concentration attributable to wildfire smoke PM2.5"),
         y = expression("Concentration (\u00b5g/m"^3*")"),
         x = "") +
    scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
    geom_hline(aes(yintercept = 0)) +
    theme_minimal() +
    scale_color_manual(values = pm_pal) +
    scale_fill_manual(values = pm_pal) +
    facet_wrap(~factor(species_long, levels = rev(limits)), scales = 'free_y', ncol = 4) +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0),
          axis.text = element_text(),
          axis.title = element_text(),
          legend.text = element_text(),
          legend.title = element_blank(),
          axis.ticks.x = element_line(linetype = 1),
          axis.ticks.length=unit(-0.2, "cm"),
          panel.grid = element_blank(),
          legend.position = "bottom",
          axis.line.y = element_line(color = "grey10"))
  avg_mon_pred_spec_plot
  
  # save file
  ggsave(
    filename = 'SIFig6_ALLsmoke_attributable_fraction_trend_regional_model.png',
    plot = avg_mon_pred_spec_plot,
    path = file.path(results_fp, 'SI Figs'),
    scale = 1,
    width = 12,
    height = 8,
    dpi = 320)
  
  # save file
  ggsave(
    filename = 'SIFig6_ALLsmoke_attributable_fraction_trend_regional_model.pdf',
    plot = avg_mon_pred_spec_plot,
    path = file.path(results_fp, 'SI Figs'),
    scale = 1,
    width = 12,
    height = 8,
    dpi = 320)
  
 # PLOT STAKCED AREA --------------------------------------------
  limits = c(
    "Organics",
    #"Organic Carbon (OC)", "Elemental Carbon (EC)",
    "Halogens",
    # "Bromine (Br)", "Chlorine (Cl)", 
    "Nonmetals",
    # "Phosphorus (P)","Sulfur (S)", "Sulfate (SO4)",  "Nitrate (NO3)", "Selenium (Se)", 
    "Other metals",
    #"Aluminum (Al)", "Lead (Pb)",
    "Metalloids",
    #"Silicon (Si)", "Arsenic (As)",
    "Transition metals",
    #"Manganese (Mn)", "Zinc (Zn)", "Titanium (Ti)", "Iron (Fe)", "Copper (Cu)", "Vanadium (V)", "Nickel (Ni)", "Chromium (Cr)",
    "Alkali metals",
    #"Potassium (K)", "Rubidium (Rb)", "Sodium (Na)",
    "Alkaline-earth metals")
    #"Strontium (Sr)","Calcium (Ca)",  "Magnesium (Mg)")
  # change to annual
  ann_plot_df <- plot_df %>%
    group_by(species_type, species_long, year, label) %>% 
    summarise(avg_ann_conc = mean(conc)) %>% 
    ungroup() %>% 
    mutate(yr = as.numeric(year)) %>% 
    mutate(species_type = factor(species_type, levels = rev(limits)))

  # Ensure the palette has exactly 27 colors (trim or repeat if needed)
  # species_colors <- c('grey20',"#899DA4", # organics
  #                     "#CAB38C", '#A67B5B', # "Halogens" 
  #                     '#62E92C', '#96BD42','#2EB173',"#0B775E",'#01401F',#  "Nonmetals"
  #                     "#46ACC8", 'steelblue', # "Other metals"
  #                     'yellow',"#E2D200", # "Metalloids"
  #                     '#FFA500' ,'#FF8C00','#FF7F50','orangered','#FF4500','#E67E22','#D2691E','#C03C21',
  #                     '#B76281', '#A12D4C', '#6D1D25',
  #                     '#DB9AFD', '#9833BD', '#872C7B')
  
avg_ann_TOC_plot <- ggplot(ann_plot_df) +
    geom_bar(aes(x = year,
                  y = avg_ann_conc,
                  fill = species_type,
                  color = species_type),
             stat = "identity", 
             position = "stack", alpha = 1) +
    labs(
      title = expression("Chemical species concentration attributable to wildfire smoke PM2.5"),
      y = expression("Concentration (\u00b5g/m"^3*")"), x = "") +
    theme_minimal() +
    facet_wrap(~label, scale = 'free_y') + 
    scale_fill_manual(values = spec_pal) +  # Use unique colors for each species
    scale_color_manual(values = spec_pal) + # Same for border colors
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0),
      axis.text = element_text(),
      axis.title = element_text(),
      legend.text = element_text(),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      axis.line.y = element_line(color = "grey10"),
      axis.line.x = element_line(color = "grey10") # Add bottom axis line
    ) +
  guides(color = 'none') +
  guides(fill = 'none')
  
avg_ann_TOC_plot
  
# save file
ggsave(
  filename = 'TOC_attributable_bar_plot_by_spec_type.pdf',
  plot = avg_ann_TOC_plot,
  path = file.path(results_fp, 'SI Figs'),
  scale = 1,
  width = 3.25,    # Set width to 3.25 inches
  height = 1.75,   # Set height to 1.75 inches
  dpi = 320
)
  
  
  # PLOT FRACTION ----------------------------------------------------------
  # prep df for plotting
  att_frac <- plot_df %>% 
    pivot_wider(names_from = 'label', values_from = 'conc') %>% 
    mutate(AF = 100*(`Attributable concentration to smoke PM2.5`/`Observed Concentration in Total PM2.5`))  %>% 
    filter(AF > 0) %>% 
    mutate(AF = ifelse(AF > 100, 100, AF))
  
  
  anno <- att_frac %>% 
    distinct(species, sig, year_trend) %>% 
    mutate(label = paste0(round(year_trend*100, digits = 1), '%')) %>% 
    mutate(label_sig = paste0(label, "-", sig)) %>% 
    filter(sig != 'cannot detect increase')
  
  df <- att_frac %>% 
    left_join(anno)
  
  datebreaks <- seq(as.Date("2006-01-01"), as.Date("2020-12-01"), by = "48 month")
  
  limits = c(
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
    "Manganese (Mn)", "Zinc (Zn)", "Titanium (Ti)", "Iron (Fe)", "Copper (Cu)", "Vanadium (V)", "Nickel (Ni)", "Chromium (Cr)",
    # "Alkali metals"
    "Potassium (K)", "Rubidium (Rb)", "Sodium (Na)",
    # "Alkaline-earth metals"
    "Strontium (Sr)","Calcium (Ca)",  "Magnesium (Mg)")
  
  # # PLOTTING
  att_frac_plot <- ggplot(df) +
    geom_density(aes(x = mon_yr,
                     y = AF,
                     fill = species_type,
                     color = species_type),
                 stat = "identity", alpha = 0.6) +
    labs(y = "% of concentration attributable to wildfire smoke",
         x = "") +
    scale_x_date(labels = date_format("%Y"), breaks = as.Date(datebreaks, format = "%Y-%m-%d")) +
    geom_hline(aes(yintercept = 0)) +
    # add trend line
    geom_smooth(aes(x = mon_yr, y = AF, group = species_long),
                method = "lm", se = FALSE, color = "black", alpha = .6, linewidth = .3) +
    geom_text(aes(x = as.Date('2015-01-01'), y = 50, label = paste(label, "/ year")),
              vjust = -0.5, hjust = 1, color = "black", size = 3, fontface = "plain") + 
    theme_minimal() +
    scale_color_manual(values = spec_pal, 
                       limits = rev(c("Organics",
                                      "Halogens",
                                      "Nonmetals",
                                      "Other metals",
                                      "Metalloids",
                                      "Transition metals",
                                      "Alkali metals",
                                      "Alkaline-earth metals")))  +
    scale_fill_manual(values = spec_pal, 
                      limits = rev(c("Organics",
                                     "Halogens",
                                     "Nonmetals",
                                     "Other metals",
                                     "Metalloids",
                                     "Transition metals",
                                     "Alkali metals",
                                     "Alkaline-earth metals"))) +
    #facet_grid(cols = vars(species_type), rows= vars(species_long)) +
    facet_wrap(~factor(species_long, levels = rev(limits)), ncol = 7) +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 2),
          axis.text = element_text(),
          axis.title = element_text(),
          legend.text = element_text(),
          legend.title = element_blank(),
          axis.ticks.x = element_line(linetype = 1),
          axis.ticks.length=unit(-0.2, "cm"),
          panel.grid = element_blank(),
          legend.position = "bottom",
          axis.line.y = element_line(color = "grey10")) +
    guides(color = 'none') +
    guides(legend = 'none') +
    guides(fill = 'none')
  att_frac_plot
  
  # save file
  ggsave(
    filename = 'Fig3_attributable_fraction_trend.png',
    plot = att_frac_plot,
    path = file.path(results_fp, 'Fig3'),
    scale = 1,
    width = 12,
    height = 8,
    dpi = 320)
  
  # save file
  ggsave(
    filename = 'Fig3_attributable_fraction_trend.pdf',
    plot = att_frac_plot,
    path = file.path(results_fp, 'Fig3'),
    scale = 1,
    width = 12,
    height = 8,
    dpi = 320)
  
  # # PLOTTING PIE CHART
  # Aggregate the data
  df_pie <- df %>%
    group_by(species) %>%
    summarise(total_AF = sum(AF, na.rm = TRUE)) %>%
    mutate(percentage = total_AF / sum(total_AF) * 100)
  
  # Create the pie chart
  att_frac_pie <- ggplot(df_pie, aes(x = "", y = total_AF, fill = species)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    labs(
      fill = "Species Type",
      title = "% of Concentration Attributable to Wildfire Smoke"
    ) +
    scale_fill_manual(values = spec_pal,
                      limits = rev(c("Organics",
                                     "Halogens",
                                     "Nonmetals",
                                     "Other metals",
                                     "Metalloids",
                                     "Transition metals",
                                     "Alkali metals",
                                     "Alkaline-earth metals"))) +
    geom_text(aes(label = round(percentage, 1)), # Label with percentages
              position = position_stack(vjust = 0.5), size = 4, color = "white") +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      legend.position = "bottom"
    )
  
  att_frac_pie
  # save file
  ggsave(
    filename = 'TOC_attributable_fraction_PIECHART.pdf',
    plot = att_frac_pie,
    path = file.path(results_fp, 'SI Figs'),
    scale = 1,
    width = 6,
    height = 6,
    dpi = 320)
  
  
  return(plot_df)
}
 