# Emma Krasovich Southworth, emmars@stanford.edu
# Last Updated: Oct 2, 2023
# Description: conduct a sensitivty analysis for our data with different objectives



# what does 0 concentration actually mean in these cases?
# mutate_if(is.numeric, ~ifelse(. == 0, NA, .)) %>% 