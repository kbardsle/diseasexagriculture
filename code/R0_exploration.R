# --------------------------------------------------------
# check to make sure R0 values for agricultural worker population within reasonable range
# 29 Apr 2024
# Katie Bardsley
# --------------------------------------------------------

# load packages
library(tidyverse)

# get path to the script directory
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# Set the working directory to the parent directory of the script directory
setwd(file.path(script_dir, ".."))

# read in data
demographic_data <- read_csv("data/migrants_merger.csv")

clean_demo_data <- demographic_data %>%
  mutate(proportion_crowded = CROWDED1.1, proportion_w_kids = (1-HHKID.0)) %>%
  select(c(State, FY, Category, Value, proportion_crowded, proportion_w_kids)) %>%
  filter(FY == 2017) %>%
  # strip commas from Value column
  mutate(Value = str_replace(Value, ",", "")) %>%
  # sum together values for two types of non-migrants
  group_by(State, FY, Category) %>%
  summarize(Value = sum(as.numeric(Value)), 
            proportion_crowded = median(as.numeric(proportion_crowded)), 
            proportion_w_kids = median(as.numeric(proportion_w_kids))) %>%
  # take weighted mean for migrant and non migrant demographic variables
  group_by(State) %>%
  summarize(proportion_crowded = weighted.mean(proportion_crowded, Value),
            proportion_w_kids = weighted.mean(proportion_w_kids, Value))

# define functions

###################################
# FUNCTION: calc_R0
# description: calculate R0 from demographic data
# inputs: prop_crowding, prop_kids, Ba, gamma
# outputs: R0
# ---------------------------------
calc_R0 <- function(prop_crowding, prop_kids, Ba, xi, eta, gamma){
  B <- Ba + xi*prop_crowding + eta*prop_kids
  R0 <- B/gamma
  return(R0)
}
###################################

clean_demo_data$R0 <- calc_R0()