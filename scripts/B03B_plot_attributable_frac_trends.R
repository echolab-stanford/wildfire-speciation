# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: April 6, 2023
# Description: make a time series plot for each species

# grid_fp = file.path(data_fp, 'intermediate/10km_grid_wgs84.shp')
# loadd(c(pred_regional_attributable_preds, spec_pal), cache = drake::drake_cache(".drake"))

plot_attributable_frac_trends <- function(pred_regional_attributable_preds, spec_pal) {

  # prep df for plotting
att_frac <- pred_regional_attributable_preds %>% 
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
  filename = 'Fig3_attributable_fraction_trend_alt.png',
  plot = att_frac_plot,
  path = file.path(results_fp, 'Fig3'),
  scale = 1,
  width = 12,
  height = 8,
  dpi = 320)

# save file
ggsave(
  filename = 'Fig3_attributable_fraction_trend_alt.pdf',
  plot = att_frac_plot,
  path = file.path(results_fp, 'Fig3'),
  scale = 1,
  width = 12,
  height = 8,
  dpi = 320)

return(df)
}
