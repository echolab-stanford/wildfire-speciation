#library
library(ggplot2)
library(dplyr)

#setup
setwd("C:/Users/ayako/OneDrive/Documents/echolab/wip-ayako/wildfire-speciation/scripts")
data_path <- "G:/Shared drives/echolab data/wildfire_speciation/intermediate/"
fig_path <- "C:/Users/ayako/OneDrive/Documents/echolab/wip-ayako/wildfire-speciation/figures"

#import data
data <- read.csv(paste(data_path, "pm_plume_speciation_at_sites.csv", sep = ''))
#file includes row names, so drop the column 
data <- data[,-1]
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


vars <- data %>% select(RC_PM2, NH4, tot_metals, ammNO3,ammSO4,
                        AL, AS, BR, CA, CR, CU, FE, K, MG, MN, NI,
                        NO3, N2, P, PB, RB, S, SE, SI, SOIL, SO4,
                        SR, V, ZN, ZR, MF) %>% names()

#log transformation
data_log <- data
data_log[vars] <- log(data_log[vars])

#all region
for(i in vars) {
  g <- ggplot(data, aes(x = get(i))) +
    geom_histogram(aes(y = ..density..),
                   binwidth = .5,
                   colour = "black", fill = "white") +
    geom_density(alpha = .2, fill = "#FF6666") +
    xlab(i)
  print(g)
  ggsave(path = paste(fig_path,"/05_hist_all_region",sep = ''), 
         filename = paste0(i, ".png"), device = "png",
         width = 5, height = 3, dpi = 300)
}


for(i in vars) {
  g <- ggplot(data_log, aes(x = get(i))) +
    geom_histogram(aes(y = ..density..),
                   binwidth = .5,
                   colour = "black", fill = "white") +
    geom_density(alpha = .2, fill = "#FF6666") +
    xlab(i)
  print(g)
  ggsave(path = paste(fig_path,"/06_hist_all_region_log",sep = ''), 
         filename = paste0(i, ".png"), device = "png",
         width = 5, height = 3, dpi = 300)
}

#each region
for(i in vars) {
g <- ggplot(data, aes(x = get(i))) +
  geom_histogram(aes(y = ..density..),
                 binwidth = .5,
                 colour = "black", fill = "white") +
  geom_density(alpha = .2, fill = "#FF6666") +
  facet_grid(abb_region ~ .) +
  xlab(i)
  print(g)
  ggsave(path = paste(fig_path,"/07_hist_by_region",sep = ''), 
       filename = paste0(i, ".png"), device = "png",
       width = 5, height = 8, dpi = 300)
}

for(i in vars) {
  g <- ggplot(data_log, aes(x = get(i))) +
    geom_histogram(aes(y = ..density..),
                   binwidth = .5,
                   colour = "black", fill = "white") +
    geom_density(alpha = .2, fill = "#FF6666") +
    facet_grid(abb_region ~ .) +
    xlab(i)
  print(g)
  ggsave(path = paste(fig_path,"/08_hist_by_region_log",sep = ''), 
         filename = paste0(i, ".png"), device = "png",
         width = 5, height = 8, dpi = 300)
}
