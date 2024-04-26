# import libraries
library(tidyverse)

# set working directory
dirname(rstudioapi::getActiveDocumentContext()$path)

# define equation to calculate beta for agricultural population
# input: base beta to modulate, list of coefficients, list of demographic variables of interest
# output: beta (numeric)
calc_beta <- function(beta_not, coeffs, vars){
  # check how many variables were passed
  num_vars <- length(vars)
  
  # multiply every variable by its associated coefficient
  partial_sums <- lapply(1:num_vars, function(i){
    curr_sum <- coeffs[[i]]*vars[[i]]
  })
  
  # add up all the partial sums and add on beta not
  return(sum(partial_sums) + beta_not)
}

# define function to determine best coefficient values
# input: desired value for r not, gamma value to use in r not calculation, base beta not value, list of variables (each item is a list of values for the variable for each location)
# output: best parameter values to achieve the desired r not based on variables passed
parameterize <- function(r_not, gamma, beta_not, vars){
  # define all values to test for each coefficient
  coeff_vals <- seq(0, 1, 0.1)
  
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
  all_r_nots <- lapply(1:nrow(all_combos), function(i){
    curr_coeffs <- list(all_combos[i,])
    # calculate r not for each location
    curr_r_nots <- lapply(all_location_vars, function(v){
      (calc_beta(beta_not, curr_coeffs, v)/gamma)
    })
    # take mean r not across all locations
    mean(curr_r_nots)
  })
  
  # calculate distance of each r not value from desired value (stress)
  stress <- abs(all_r_nots - r_not)
  
  # return combination of coefficients with lowest stress
  best_coeffs <- list(all_combos[which.min(stress),])
}



