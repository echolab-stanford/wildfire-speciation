#library
library(ggplot2)
library(dplyr)

#setup
setwd("C:/Users/ayako/Documents/StanfordPhD/EchoLab/wildfire-speciation/scripts")
data_path <- "G:/Shared drives/echolab data/wildfire_speciation/intermediate/"
fig_path <- "C:/Users/ayako/Documents/StanfordPhD/EchoLab/wildfire-speciation/figures"

#import data
data <- read.csv(paste(data_path, "pm_plume_speciation_at_sites_tempAK.csv", sep = ''))
#file includes row names, so drop the column 
#data <- data[,-1]
data$smoke_day <- as.factor(data$smoke_day)
#add a short region name
data <- data %>% 
  mutate(abb_region = case_when(
    region == "rocky_mountain" ~ "RM",
    region =="midwest" ~ "MW",
    region =="southwest" ~ "SW",
    region =="noncontiguous" ~ "Non",
    region =="pacific" ~ "PAC",
    region =="northeast" ~ "NE",
    region =="southeast" ~ "SE"
))


summary(data)
colnames(data)
vars <- data %>% select(MF, RC_PM2, NH4, tot_metals, ammNO3, ammSO4, AL, AS, BR, 
                        CA, EC, OC, CHL, CL, CR, CU, FE, K, MG, MN, MO, NA., NI,
                        NO3, N2, P, PB, RB, S, SE, SI, SOIL, SO4, SR, V, ZN, ZR
                        )%>% names()

#log transformation
data_log <- data
data_log[vars] <- log(data_log[vars])

#all region
for(i in vars) {
  g <- ggplot(data, aes(x = smoke_day, y = get(i))) +
    geom_boxplot() +
    ylab(i) +
    xlab("Smoke Day")
  print(g)
  ggsave(path = paste(fig_path,"/01_boxplots_all_region",sep = ''), 
         filename = paste0(i, ".png"), device = "png",
         width = 5, height = 3, dpi = 300)
}

for(i in vars) {
  g <- ggplot(data_log, aes(x = smoke_day, y = get(i))) +
    geom_boxplot() +
    ylab(i) +
    xlab("Smoke Day")
  print(g)
  ggsave(path = paste(fig_path,"/02_boxplot_all_region_log",sep = ''), 
         filename = paste0(i, ".png"), device = "png",
         width = 5, height = 3, dpi = 300)
}


#each region
for(i in vars) {
  g <- ggplot(data, aes(x = abb_region, y = get(i), fill = smoke_day)) +
    geom_boxplot() +
    ylab(i) 
  print(g)
  ggsave(path = paste(fig_path,"/03_boxplot_by_region",sep = ''), 
         filename = paste0(i, ".png"), device = "png",
         width = 5, height = 3, dpi = 300)
}

for(i in vars) {
  g <- ggplot(data_log, aes(x = abb_region, y = get(i), fill = smoke_day)) +
    geom_boxplot() +
    ylab(i) 
  print(g)
  ggsave(path = paste(fig_path,"/04_boxplot_by_region_log",sep = ''), 
         filename = paste0(i, ".png"), device = "png",
         width = 5, height = 3, dpi = 300)
}

#CSN vs IMPROVE
for(i in vars) {
  g <- ggplot(data_log, aes(x = abb_region, y = get(i), fill = Dataset)) +
    geom_boxplot() +
    ylab(i) 
  print(g)
  ggsave(path = paste(fig_path,"/12_boxplot_by_monitor_log",sep = ''), 
         filename = paste0(i, ".png"), device = "png",
         width = 5, height = 3, dpi = 300)
}
