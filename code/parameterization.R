# parameterize model with demographic data from the general population

# import libraries
library(tidyverse)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

# define equation to calculate beta for agricultural population
# input: base beta to modulate, list of coefficients, list of demographic variables of interest
# output: beta (numeric)
calc_beta <- function(beta_not, coeffs, vars){
  # check how many variables were passed
  num_vars <- length(vars)
  # print("num vars")
  # print(num_vars)
  # print("coeffs")
  # print(coeffs)
  
  # multiply every variable by its associated coefficient
  partial_sums <- do.call(rbind, lapply(1:num_vars, function(i){
    # print("curr index")
    # print(i)
    curr_sum <- coeffs[[i]]*vars[[i]]
  }))
  
  # add up all the partial sums and add on beta not
  return(sum(partial_sums) + beta_not)
}

# define function to determine best coefficient values
# input: desired value for r not, gamma value to use in r not calculation, base beta not value, list of variables (each item is a list of values for the variable for each location)
# output: best parameter values to achieve the desired r not based on variables passed
parameterize <- function(r_not, gamma, beta_not, vars){
  # define all values to test for each coefficient
  coeff_vals <- seq(0, 1, 0.01)
  
  # check how many variables there are
  num_vars <- length(vars)
  # construct dataframe of every possible combination of coefficient values
  all_combos <- expand.grid(replicate(num_vars, coeff_vals, simplify = FALSE))
  
  # check how many locations there are
  num_locations <- length(vars[[1]])
  
  # pull together all variables for each location
  all_locations_vars <- lapply(1:num_locations, function(l){
    # gather all variables for current location
    curr_location_vars <- lapply(1:num_vars, function(v){
      curr_var <- vars[[v]][[l]]
    })
  })
  
  # calculate r not for every combination of coefficients
  all_r_nots <- do.call(rbind, lapply(1:nrow(all_combos), function(i){
    curr_coeffs <- list(all_combos[i,])
    # calculate r not for each location
    curr_r_nots <- do.call(rbind, lapply(all_locations_vars, function(v){
      (calc_beta(beta_not, curr_coeffs[[1]], v)/gamma)
    }))
    # take mean r not across all locations
    mean(curr_r_nots)
  }))
  
  # calculate distance of each r not value from desired value (stress)
  stress <- abs(all_r_nots - r_not)
  
  print(min(stress))
  
  # return combination of coefficients with lowest stress
  best_coeffs <- list(all_combos[which.min(stress),])
}

# read in demographic data
demographic_data_raw <- read_csv("data/general_population_demographics.csv")

# normalize demographic data
demographic_data_normalized <- demographic_data_raw %>%
  mutate(proportion_crowded = percent_crowded/100,
         proportion_with_children = percent_with_children/100)

demographic_vars <- list(demographic_data_normalized$proportion_crowded, demographic_data_normalized$proportion_with_children)

# run parameterization code
coefficient_vals <- parameterize(r_not = 1.6,
                                 gamma = 0.125,
                                 beta_not = 0.1,
                                 vars = demographic_vars)


