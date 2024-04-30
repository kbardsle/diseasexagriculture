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
data_2017 <- read_csv("data/2017_pop_demo_data_agricultural.csv")

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

data_2017$R0 <- calc_R0(data_2017$proportion_crowded, data_2017$proportion_w_kids, 0.1, 0.92, 0.25, 0.125)

summary(data_2017$R0)
weighted.mean(data_2017$R0, data_2017$POP3)

ggplot(data=data_2017, mapping=aes(x=R0)) + geom_histogram(bins=50) + theme_bw()
