# Ayako Kawano, Emma Krasovich Southworth
# Description: This 


#import data
data <- read.csv(paste(data_path, "pm_plume_speciation_at_sites.csv", sep = ''))

# Group species
columns_list <-
  c("V", "PB", "CR",
    "FE", "CU", "MN", "ZN", "NI",
    "EC", "OC",
    "S")

df <- data %>% select("month", "totPM2.5", "smokePM", "V", "PB", "CR",
                      "FE", "CU", "MN", "ZN", "NI","EC", "OC","S") %>% 
  rename("Smoke PM2.5" = "smokePM", "Total PM2.5" = "totPM2.5")


monthly <- df %>% group_by(month) %>% 
  summarise(across(everything(), mean, na.rm = T)) %>% 
  reshape2::melt(id = "month") %>% 
  mutate(
        Category = ifelse(variable=="V" | variable=="PB" | variable=="CR", "Toxic metal",
                  ifelse(variable=="FE" | variable=="CU" |
                           variable=="MN" | variable=="ZN" | variable=="NI", "Transition metal",
                         ifelse(variable=="EC" | variable=="OC", "Organic",
                                ifelse(variable=="S", "Toxicity potentiator", 
                                       ifelse(variable=="Total PM2.5" | variable=="Smoke PM2.5", "PM2.5",NA)))))
  )

#monthly <- monthly %>% add_row(month = 13, variable = NA,  value= NA, Category = NA)

#df including only species, not PM
species <- monthly %>% filter(Category != "PM2.5")

#plot pct change setting January as a baseline----
pct <- monthly %>% filter(variable != "Smoke PM2.5") %>% 
  mutate(pct_change = NA)
pct <- pct %>% add_row(month = 13, variable = NA,  value= NA, Category = NA, pct_change = NA)


pct <- pct %>%
  group_by(variable) %>% 
  mutate(pct_change = (value/lag(value) - 1) * 100)

pct <- pct %>% mutate(pct_change = replace(pct_change,is.na(pct_change),0))
pct[145, "pct_change"] <- NA

#pct <- pct %>% mutate(variable = fct_relevel(variable,
#                                             "Total PM2.5", "OC", "EC",
#                                             "V", "PB", "CR",
#                                             "FE", "CU", "MN", "ZN", "NI","S"))

pct <- pct %>% mutate(Category = fct_relevel(Category,
                                             c("PM2.5", "Organic", "Toxic metal",
                                             "Transition metal", "Toxicity potentiator")))
levels(pct$Category)
cbp1 <- c("firebrick3","goldenrod1", "goldenrod3","goldenrod4",
          "lightcyan2", "gray80", "lightblue3", "darkslategray3","gray40",
          "palegreen3",
          "palegreen4", "deepskyblue1")
data_ends <- pct %>% filter(month == 12)

ggplot(pct, aes(x = month, y = pct_change, group = variable)) +
  geom_line(size = 0.8, aes(color = variable)) +
  ylab("Change (%)") +
  xlab("") +
  ggtitle("Seasonal differences in concentrations of chemical species") +
  scale_x_continuous(breaks = seq(1,12,1),
                     labels = c("Jan", "Feb", "Mar", "Apr",
                                "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  scale_y_continuous(limits = c(-50,100) , breaks=c(-50, -25, 0, 25, 50, 75, 100),
                     labels=c("-50", "-25", "0","25", "50", "75", "100")) +
  #geom_point(size = 1.2, aes(color = Category)) +
  geom_point(data = data_ends, size =1.5,
             aes(x = month, y = pct_change, color = variable), show.legend = FALSE)+
  #geom_text(data = data_ends,
  #          aes(x =12, 
  #              y = pct_change,  
  #              label = variable), nudge_x = 0.2, show.legend = FALSE, hjust = 0, size = 3)+
  #scale_linetype_manual(values=c("solid", "longdash", "dotted","twodash", "dashed")) +
  theme_minimal() +
  #coord_cartesian() +
  scale_color_manual(values = cbp1, na.translate = F) +
  theme(plot.title = element_text(size = 12, family ="serif", face = "bold", hjust = 0.4),
        axis.text = element_text(size = 12, family ="serif"),
        axis.title = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(vjust = 0),
        legend.text = element_text(size = 10, family = "serif"),
        legend.title = element_blank(),
        axis.ticks.x = element_line(linetype = 1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.y = element_line(color = "gray90")) 
  #axis.line.x = element_line(size = 0.5, colour = "black", linetype=1)) +
  #geom_hline(aes(yintercept = 0))


#plot all region, all years, monthly means for each species including PM2.5 ----
cbp1 <- c("firebrick3", "#BDB76B", "orange1", "deepskyblue1", "gray80")

data_ends <- monthly %>% filter(month == 12)

pdf(file="G:/Shared drives/echolab data/wildfire_speciation/figures/seasonality/seasonality_all.pdf")  
ggplot(monthly, aes(x = month, y = value, group = variable)) +
  geom_line(size = 0.8, aes(color = Category)) +
  ylab("Concentration (ug/m3)") +
  xlab("") +
  ggtitle("Seasonal differences in concentrations of chemical species") +
  scale_x_continuous(breaks = seq(1,13,1),
                     labels = c("Jan", "Feb", "Mar", "Apr",
                     "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "")) +
  scale_y_continuous(trans='log10', limits = c(10^-5,9) , breaks=c(10^-5, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 4.0, 8.0),
                     labels=c("0","0.0005", "0.001", "0.005", "0.01", "0.05", "0.1", "0.5", "1.0", "2.0", "4.0", "8.0")) +
  #geom_point(size = 1.2, aes(color = Category)) +
  geom_point(data = data_ends, size =1.5,
             aes(x = month, y = value, color = Category), show.legend = FALSE)+
  geom_text(data = data_ends,
            aes(x =12, 
                y = value,  
                label = variable), nudge_x = 0.2, show.legend = FALSE, hjust = 0, size = 3)+
  #scale_linetype_manual(values=c("solid", "longdash", "dotted","twodash", "dashed")) +
  theme_minimal() +
  coord_cartesian() +
  scale_color_manual(values = cbp1, na.translate = F) +
  theme(plot.title = element_text(size = 12, family ="serif", face = "bold", hjust = 0.4),
        axis.text = element_text(size = 12, family ="serif"),
        axis.title = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(vjust = 0),
        legend.text = element_text(size = 10, family = "serif"),
        legend.title = element_text(size = 10, family = "serif"),
        axis.ticks.x = element_line(linetype = 1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.y = element_line(color = "gray90")) +
        #axis.line.x = element_line(size = 0.5, colour = "black", linetype=1)) +
  geom_hline(aes(yintercept = 10^-5))
dev.off()



#plot all region, all years, monthly means for each species without PM2.5 ----
cbp2 <- c("firebrick3", "orange1", "deepskyblue1", "gray80")

data_ends <- species %>% filter(month == 12)


ggplot(species, aes(x = month, y = value, group = variable)) +
  geom_line(size = 0.8, aes(color = Category)) +
  ylab("Concentration (ug/m3)") +
  xlab("") +
  ggtitle("Seasonal differences in concentrations of chemical species") +
  scale_x_continuous(breaks = seq(1,13,1),
                     labels = c("Jan", "Feb", "Mar", "Apr",
                                "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "")) +
  scale_y_continuous(trans='log10', limits = c(10^-5,9) , breaks=c(10^-5, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0),
                     labels=c("0","0.0005", "0.001", "0.005", "0.01", "0.05", "0.1", "0.5", "1.0", "2.0")) +
  #geom_point(size = 1.2, aes(color = Category)) +
  geom_point(data = data_ends, size =1.5,
             aes(x = month, y = value, color = Category), show.legend = FALSE)+
  geom_text(data = data_ends,
            aes(x =12, 
                y = value,  
                label = variable), nudge_x = 0.2, show.legend = FALSE, hjust = 0, size = 3)+
  #scale_linetype_manual(values=c("solid", "longdash", "dotted","twodash", "dashed")) +
  theme_minimal() +
  coord_cartesian() +
  scale_color_manual(values = cbp2, na.translate = F) +
  theme(plot.title = element_text(size = 12, family ="serif", face = "bold", hjust = 0.4),
        axis.text = element_text(size = 12, family ="serif"),
        axis.title = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(vjust = 0),
        legend.text = element_text(size = 10, family = "serif"),
        legend.title = element_text(size = 10, family = "serif"),
        axis.ticks.x = element_line(linetype = 1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.y = element_line(color = "gray90")) +
  #axis.line.x = element_line(size = 0.5, colour = "black", linetype=1)) +
  geom_hline(aes(yintercept = 10^-5))


#density plot version: plot all region, all years, monthly means for each species without PM2.5 ----
#cbp2 <- c("firebrick3", "orange1", "deepskyblue1", "gray80")

data_ends <- species %>% filter(month == 12)

ggplot(species, aes(x = month, y = value, fill = variable, color = fct_rev(Category))) +
  geom_density(stat = "identity", position = "stack", alpha = 0.6) +  
  ylab("Concentration (ug/m3)") +
  xlab("") +
  ggtitle("Seasonal differences in concentrations of chemical species") +
  scale_x_continuous(breaks = seq(1,13,1),
                     labels = c("Jan", "Feb", "Mar", "Apr",
                                "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "")) +
  scale_y_continuous(trans='log10', limits = c(10^-5,9) , breaks=c(10^-5, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 4.0),
                     labels=c("0","0.0005", "0.001", "0.005", "0.01", "0.05", "0.1", "0.5", "1.0", "2.0", "4.0")) +
  #geom_point(size = 1.2, aes(color = Category)) +
  geom_point(data = data_ends, size =1.5,
             aes(x = month, y = value, color = Category), show.legend = FALSE)+
  geom_text(data = data_ends,
            aes(x =12, 
                y = value,  
                label = variable), nudge_x = 0.2, show.legend = FALSE, hjust = 0, size = 3)+
  #scale_linetype_manual(values=c("solid", "longdash", "dotted","twodash", "dashed")) +
  theme_minimal() +
  coord_cartesian() +
  scale_color_manual(values = cbp2, na.translate = F) +
  theme(plot.title = element_text(size = 12, family ="serif", face = "bold", hjust = 0.4),
        axis.text = element_text(size = 12, family ="serif"),
        axis.title = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(vjust = 0),
        legend.text = element_text(size = 10, family = "serif"),
        legend.title = element_text(size = 10, family = "serif"),
        axis.ticks.x = element_line(linetype = 1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.y = element_line(color = "gray90")) +
  #axis.line.x = element_line(size = 0.5, colour = "black", linetype=1)) +
  geom_hline(aes(yintercept = 10^-5))



#+ geom_text_repel(
#  aes(label = variable), data = data_ends, size = 2.5, 
#  nudge_x = 0.28, nudge_y = 0,direction = "y", hjust = "right", family = "serif"
#) 

#plot all region, all years, monthly means for PM2.5 ----
cbp3 <- c("gray80", "firebrick3")

pm25 <- monthly %>% filter (Category == "PM2.5")
data_ends <- pm25 %>% filter(month == 12)

pm25$variable <- as.character(pm25$variable, levels = c("Total PM2.5", "Smoke PM2.5"))

ggplot(pm25) +
  geom_density(aes(x = month, y = value, fill = fct_rev(variable), color = fct_rev(variable)), stat = "identity", position = "stack", alpha = 0.6) +  
  ylab("Concentration (ug/m3)") +
  xlab("") +
  ggtitle("Seasonal differences in concentrations of PM2.5") +
  scale_x_continuous(breaks = seq(1,13,1),
                     labels = c("Jan", "Feb", "Mar", "Apr",
                                "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "")) +
  scale_y_continuous(limits = c(0,10) , breaks=c(0, 2, 4, 6, 8, 10),
                     labels=c("0","2.0", "4.0", "6.0", "8.0", "10.0")) +
  #geom_point(size = 1.2, aes(color = Category)) +
  #geom_point(data = data_ends, size =1.5,
   #          aes(x = month, y = value, color = variable), show.legend = FALSE)+
  #geom_text(data = data_ends,
  #          aes(x =12, 
  #              y = value,  
  #              label = variable), nudge_x = 0.2, show.legend = FALSE, hjust = 0, size = 3)+
  #scale_linetype_manual(values=c("solid", "longdash", "dotted","twodash", "dashed")) +
  theme_minimal() +
  #coord_cartesian() +
  scale_color_manual(values = cbp3, na.translate = F) +
  scale_fill_manual(values = cbp3, na.translate = F) +
  theme(plot.title = element_text(size = 12, family ="serif", face = "bold", hjust = 0.4),
        axis.text = element_text(size = 12, family ="serif"),
        axis.title = element_text(size = 12, family = "serif"),
        axis.text.x = element_text(vjust = 1),
        legend.text = element_text(size = 10, family = "serif"),
        legend.title = element_blank(),
        axis.ticks.x = element_line(linetype = 1),
        axis.ticks.length=unit(-0.1, "cm"),
        panel.grid = element_blank(),
        #panel.border = element_blank(),
        panel.grid.major.y = element_line(color = "gray90")) +
  #axis.line.x = element_line(size = 0.5, colour = "black", linetype=1)) +
  geom_hline(aes(yintercept = 0))






#all region, time series of monthly means from 2006-2021
#for(i in vars) {
  
#  g <- ggplot(monthly_gb, aes(x = month_yr, y = get(i))) +
#    geom_point(na.rm=TRUE, color="darkblue", size=1) +
#    (scale_x_date(breaks=date_breaks("12 months"),
#                  labels=date_format("%y"))) +
#    ylab(i) +
#    xlab("Year") +
#    ggtitle("Monthly means, 2006-2021")
#  print(g)
#  ggsave(path = paste(fig_path,"/09_monthly_means",sep = ''), 
#         filename = paste0(i, ".png"), device = "png",
#         width = 5, height = 3, dpi = 300)
#}

#all region, seasonality plot of monthly means from 2006-2021
#ggplot(monthly_gb, aes(x = month, y = RC_PM2, group = year, color = year)) +
#  geom_line()

#all region, time series of monthly means from 2006-2021
#for(i in vars) {
  
#  g <- ggplot(monthly_gb, aes(x = month, y = get(i), group = year, color = year)) +
#    geom_line() +
#    ylab(i) +
#    xlab("Month") +
#    ggtitle("Monthly means, 2006-2021")
#  print(g)
#  ggsave(path = paste(fig_path,"/10_monthly_means_seasonality",sep = ''), 
#         filename = paste0(i, ".png"), device = "png",
#         width = 5, height = 3, dpi = 300)
#}

