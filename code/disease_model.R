# Defining spatially and temporally explicit disease model using Odin
# Group rotation project - Spring 2024

# load libraries
library(odin)
library(tidyverse)

# build spatial model
sir <- odin::odin({

  # Logic for SIR Calculations:
    # if a location is newly infected (infected[i] > 0 and I[i] == 0), 
        # increase the proportion of infected by a set seed (here 0.01) to seed the infection
    # if a location was infected in an earlier time step, calculate derivatives as normal
    # if a location is still not infected, set derivatives to 0 (don't change)

  # SUSCEPTIBLE:
  deriv(S[]) <- if (infected[i] > 0 && I[i] == 0) -beta*S[i]*(I[i]+0.01) else if (infected[i] > 0 && I[i] > 0) -beta*S[i]*I[i] else 0
  # deriv(S[]) <- if (infected[i] > 0 && I[i] == 0) -beta*S[i]*(1/N[i]) else if (infected[i] > 0 && I[i] > 0) -beta*S[i]*I[i] else 0  # here edited the seed so 1 over the population size (the proportion equal to 1 infected individual)
  
  # INFECTED:
  deriv(I[]) <- if (infected[i] > 0 && I[i] == 0) beta*S[i]*(I[i]+0.01) - gamma*(I[i]+0.01) else if (infected[i] > 0 && I[i] > 0 ) beta*S[i]*I[i] - gamma*I[i] else 0
  # deriv(I[]) <- if (infected[i] > 0 && I[i] == 0) beta*S[i]*(1/N[i]) - gamma*(1/N[i]) else if (infected[i] > 0 && I[i] > 0 ) beta*S[i]*I[i] - gamma*I[i] else 0  # here edited the seed so 1 over the population size (the proportion equal to 1 infected individual)
  
  
  # RECOVERED:
  deriv(R[]) <- if (infected[i] > 0 && I[i] == 0) gamma*(I[i]+0.01) else if (infected[i] > 0 && I[i] > 0 ) gamma*I[i] else 0
  # deriv(R[]) <- if (infected[i] > 0 && I[i] == 0) gamma*(1/N[i]) else if (infected[i] > 0 && I[i] > 0 ) gamma*I[i] else 0  # here edited the seed so 1 over the population size (the proportion equal to 1 infected individual)
  
  # DEBUGGING:
  deriv(x[]) <- infected[i]
  
  # DEFINE TERMS:
  
  k_dij[,] <- exp(-dij[i,j]/ro)  # exponential kernel
  
  # matrix of the values for the equation numerator for all combinations of i and j:
  NK_num[,] <- (N[j] * if (I[j] > 0) 1 else 0)^nu * k_dij[i,j]
      # value is 0 if location not infected, in order for this value to not be included in the sum
      # because only summing up values for infected locations
  # matrix of the values for the equation denominator for all combinations of i and j:
  NK_denom[,] <- N[j]^nu * k_dij[i,j]  
  # both the numerator and denominator matrices will be further filtered before summations are taken (see below)
  
  # calculate lambda (the force of infection) (should be between 0 and infinity):
      # for the numerator, a value is included in the sum if it is row i (and not the location is infected)
      # for the denominator, a value is included in the sum if it is in row i, 
      # as long as it isn't on the diagonal (i.e., i != j)
  lambda[] <- if (I[i] == 0) beta_not + (beta_d + (school * beta_ds)) * N[i]^mu * (sum(NK_num[i,]))/((sum(NK_denom[i,]) - NK_denom[i,i])^epsilon) else 0
  
  # use lambda to calculate the probability of infection starting in a location:
  infection_prob[] <- (1-exp(-lambda[i]))
  
  # create a vector tracking whether a location has been infected
      # 1 for infected, 0 if not infected yet
      # will be marked infected if the infection probability is above a user-specified threshold
  infected[] <- if ((infection_prob[i] > infection_threshold) || (I[i] > 0)) 1 else 0

  # PRINT STATEMENTS TO HELP WITH DEBUGGING:
  
  # print("infected: {infected[1]}, {infected[2]}, {infected[3]}, {infected[4]}") #, when = infected[1] == 0)
  # print("infection probs: {infection_prob[1]}, {infection_prob[2]}, {infection_prob[3]}, {infection_prob[4]}")
  print("lambda: {lambda[1]}, {lambda[2]}, {lambda[3]}, {lambda[4]}")
  # print("NK_denom: {NK_denom[1]}, {NK_denom[2]}, {NK_denom[3]}, {NK_denom[4]}")
  print("(sum(NK_denom[1,]) - NK_denom[1,1])^epsilon: sum({NK_denom[1,]}) - {NK_denom[1,1]})^{epsilon}")

  
  
  # DEFINE USER INPUTS
  
  N[] <- user()  # vector of population sizes for locations
  dij[,] <- user()  # matrix of values for distances between i and j locations
  init_S[] <- user()  # vector of starting susceptible population proportions for each location
  init_I[] <- user()  # vector of starting infected population proportions for each location
  init_R[] <- user()  # vector of starting recovered population proportions for each location
  
  # parameters
  beta_not <- user()
  beta_d <- user()
  beta_ds <- user()
  school <- user()
  mu <- user()
  nu <- user()
  epsilon <- user()
  n_locations <- user()
  infection_threshold <- user()
  ro <- user()
  beta <- user()
  gamma <- user()
  
  # DEFINE INITIAL STATE
  initial(S[]) <- init_S[i]
  initial(I[]) <- init_I[i]
  initial(R[]) <- init_R[i]
  
  initial(x[]) <- 0  # added to help with debugging
  
  # ASSIGN DIMENSIONS
  dim(N) <- n_locations
  dim(k_dij) <- c(n_locations, n_locations)
  dim(infected) <- n_locations
  dim(dij) <- c(n_locations, n_locations)
  dim(S) <- n_locations
  dim(I) <- n_locations
  dim(R) <- n_locations
  dim(NK_num) <- c(n_locations, n_locations)
  dim(NK_denom) <- c(n_locations, n_locations)
  dim(lambda) <- n_locations
  dim(infection_prob) <- n_locations
  dim(init_S) <- n_locations
  dim(init_I) <- n_locations
  dim(init_R) <- n_locations
  
  dim(x) <- n_locations  # added to help with debugging
  
}, debug_enable=TRUE)


# TEST OUR MODEL

# build test data
populations <- c(100, 100, 100, 100)
distances_vec <- c(0, 5, 10, 8,
                   5, 0, 20, 7,
                   10, 20, 0, 15,
                   8, 7, 15, 0)*1000
distances <- matrix(data = distances_vec, nr = 4, nc = 4)

# define parameters
beta_not <- 0.0004
beta_d <- 0.77
beta_ds <- 0
school <- 0
mu <- 0.23
nu <- 0
epsilon <- 1
n_locations <- 4
infection_threshold <- 0.5   # changed this from 0.001
ro <- 96
beta <- 4
gamma <- 0  # this was 3.5, changed to 0 so nobody would recover to help with debugging

# user inputs
init_S <- c(1, 1, .99, 1)
init_I <- c(0, 0, .01, 0)
init_R <- c(0, 0, 0, 0)

t <- seq(from=0, to=10, by=1)  # change to from 100 to see if just taking a long time

model <- sir$new(beta_not=beta_not,
                 beta_d=beta_d,
                 beta_ds=beta_ds,
                 school=school,
                 mu=mu, nu=nu,
                 epsilon=epsilon,
                 n_locations=n_locations,
                 infection_threshold=infection_threshold,
                 ro=ro,
                 beta=beta,
                 gamma=gamma,
                 init_S=init_S,
                 init_I=init_I,
                 init_R=init_R,
                 N=populations,
                 dij=distances)

sol <- model$run(t)


# check if values add to 1
one_check <- sol %>%
  as.data.frame() %>%
  rename('S1' = 'S[1]',
         'I1' = 'I[1]',
         'R1' = 'R[1]') %>%
  mutate(one = S1 + I1 + R1)

one_check$one




