# Defining spatially and temporally explicit disease model using Odin
# Group rotation project - Spring 2024

# load libraries
library(odin)
library(tidyverse)
library(codep)  # for calculating great circle distances


# set working directory



# build spatial model
sir <- odin::odin({

  # Logic for SIR Calculations:
    # if a location is newly infected (infected[i] > 0 and I[i] == 0), 
        # increase the proportion of infected by a set seed (here 0.01) to seed the infection
    # if a location was infected in an earlier time step, calculate derivatives as normal
    # if a location is still not infected, set value equal to previous time step

  # Community transmission:
  
  # SUSCEPTIBLE:
  update(Sc[]) <- if (infected[i] > 0 && Ic[i] == 0) Sc[i] - beta_c*Sc[i]*(Ic[i]+seed) else if (infected[i] > 0 && Ic[i] > 0) Sc[i] - beta_c*Sc[i]*Ic[i] else Sc[i]

  # INFECTED:
  update(Ic[]) <- if (infected[i] > 0 && Ic[i] == 0) Ic[i] + beta_c*Sc[i]*(Ic[i]+seed) - gamma*(Ic[i]+seed) else if (infected[i] > 0 && Ic[i] > 0 ) Ic[i] + beta_c*Sc[i]*Ic[i] - gamma*Ic[i] else Ic[i]
  
  # RECOVERED:
  update(Rc[]) <- if (infected[i] > 0 && Ic[i] == 0) Rc[i] + gamma*(Ic[i]+seed) else if (infected[i] > 0 && Ic[i] > 0 ) Rc[i] + gamma*Ic[i] else Rc[i]

  # DEBUGGING:
  # update(x[]) <- infected[i]
  
  # Agricultural workforce transmission:
  
  # SUSCEPTIBLE:
  update(Sa[]) <- Sa[i] - ((beta_c*Ic[i])+((beta_a+xi*avg_house[i]+eta*avg_age[i])*Ia[i]))*Sa[i]
  
  # INFECTED:
  update(Ia[]) <- Ia[i] + ((beta_c*Ic[i])+((beta_a+xi*avg_house[i]+eta*avg_age[i])*Ia[i]))*Sa[i] - gamma*Ia[i]
  
  # RECOVERED:
  update(Ra[]) <- Ra[i] + gamma*Ia[i]
  
  # DEFINE TERMS:
  
  k_dij[,] <- exp(-dij[i,j]/ro)  # exponential kernel
  
  # matrix of the values for the equation numerator for all combinations of i and j:
  NK_num[,] <- (N[j] * if (Ic[j] > 0) 1 else 0)^nu * k_dij[i,j]    # now that I is no longer a derivative, does this logic still hold? (same with lambda below) **********
      # value is 0 if location not infected, in order for this value to not be included in the sum
      # because only summing up values for infected locations
  # matrix of the values for the equation denominator for all combinations of i and j:
  NK_denom[,] <- N[j]^nu * k_dij[i,j]  # for testing, nu=0, so this should just be the distance kernel
  # both the numerator and denominator matrices will be further filtered before summations are taken (see below)
  
  # calculate lambda (the force of infection) (should be between 0 and infinity):
      # for the numerator, a value is included in the sum if it is row i (and not the location is infected)
      # for the denominator, a value is included in the sum if it is in row i, 
      # as long as it isn't on the diagonal (i.e., i != j)
  NK_denom_sum[] <- sum(NK_denom[i,])  # added for debugging
  denom[] <- (sum(NK_denom[i,]) - NK_denom[i,i])^epsilon
  lambda[] <- if (Ic[i] == 0) beta_not + (beta_d + (school * beta_ds)) * N[i]^mu * (sum(NK_num[i,]))/(denom[i]) else 0
  # lambda[] <- if (I[i] == 0) beta_not + (beta_d + (school * beta_ds)) * N[i]^mu * (sum(NK_num[i,]))/((sum(NK_denom[i,]) - NK_denom[i,i])^epsilon) else 0
  
  # use lambda to calculate the probability of infection starting in a location:
  infection_prob[] <- (1-exp(-lambda[i]))
  
  # create a vector tracking whether a location has been infected
      # 1 for infected, 0 if not infected yet
      # will be marked infected if the infection probability is above a user-specified threshold
  infection_draws[] <- runif(0,1)
  infected[] <- if (infection_draws[i] <= infection_prob[i]) 1 else 0

  # PRINT STATEMENTS TO HELP WITH DEBUGGING:
  
  # print("infected: {infected[1]}, {infected[2]}, {infected[3]}, {infected[4]}") #, when = infected[1] == 0)
  print("infection probs: {infection_prob[1]}, {infection_prob[2]}, {infection_prob[3]}, {infection_prob[4]}")
  # print("lambda: {lambda[1]}, {lambda[2]}, {lambda[3]}, {lambda[4]}")
  # print("dij: {dij[1,1]}, {dij[1,2]}, {dij[1,3]}, {dij[1,4]}")
  # print("k_dij: {k_dij[1,1]}, {k_dij[1,2]}, {k_dij[1,3]}, {k_dij[1,4]}")
  # print("NK_denom: {NK_denom[1,1]}, {NK_denom[1,2]}, {NK_denom[1,3]}, {NK_denom[1,4]}")
  # print("NK_denom_sum: {NK_denom_sum[1]}, {NK_denom_sum[2]}, {NK_denom_sum[3]}, {NK_denom_sum[4]}")
  # print("denom: {denom[1]}, {denom[2]}, {denom[3]}, {denom[4]}")
  

  
  
  # DEFINE USER INPUTS
  
  N[] <- user()  # vector of population sizes for locations
  dij[,] <- user()  # matrix of values for distances between i and j locations
  # community level:
  init_Sc[] <- user()  # vector of starting susceptible population proportions for each location
  init_Ic[] <- user()  # vector of starting infected population proportions for each location
  init_Rc[] <- user()  # vector of starting recovered population proportions for each location
  # agricultural population:
  init_Sa[] <- user()  # vector of starting susceptible population proportions for each location
  init_Ia[] <- user()  # vector of starting infected population proportions for each location
  init_Ra[] <- user()  # vector of starting recovered population proportions for each location
  
  # parameters
  beta_not <- user()
  beta_d <- user()
  beta_ds <- user()
  school <- user()
  mu <- user()
  nu <- user()
  epsilon <- user()
  n_locations <- user()
  ro <- user()
  beta_c <- user()
  beta_a <- user()
  gamma <- user()
  xi <- user()
  eta <- user()
  seed <- user()
  
  # agricultural workforce community parameters
  avg_age[] <- user()
  avg_house[] <- user()
  
  # DEFINE INITIAL STATE
  initial(Sc[]) <- init_Sc[i]
  initial(Ic[]) <- init_Ic[i]
  initial(Rc[]) <- init_Rc[i]
  
  initial(Sa[]) <- init_Sa[i]
  initial(Ia[]) <- init_Ia[i]
  initial(Ra[]) <- init_Ra[i]
  
  # initial(x[]) <- 0  # added to help with debugging
  
  # ASSIGN DIMENSIONS
  dim(N) <- n_locations
  dim(k_dij) <- c(n_locations, n_locations)
  dim(infected) <- n_locations
  dim(dij) <- c(n_locations, n_locations)
  dim(Sc) <- n_locations
  dim(Ic) <- n_locations
  dim(Rc) <- n_locations
  dim(Sa) <- n_locations
  dim(Ia) <- n_locations
  dim(Ra) <- n_locations
  dim(NK_num) <- c(n_locations, n_locations)
  dim(NK_denom) <- c(n_locations, n_locations)
  dim(NK_denom_sum) <- n_locations  # added for debugging
  dim(denom) <- n_locations  # added for debugging
  dim(lambda) <- n_locations
  dim(infection_prob) <- n_locations
  dim(init_Sc) <- n_locations
  dim(init_Ic) <- n_locations
  dim(init_Rc) <- n_locations
  dim(init_Sa) <- n_locations
  dim(init_Ia) <- n_locations
  dim(init_Ra) <- n_locations
  dim(infection_draws) <- n_locations
  dim(avg_house) <- n_locations
  dim(avg_age) <- n_locations
  
  # dim(x) <- n_locations  # added to help with debugging
  
}, debug_enable=TRUE)


# TEST OUR MODEL

# build test data
populations <- c(100, 100, 100, 100)
distances_vec <- c(0, 5, 10, 80,
                   5, 0, 200, 7,
                   10, 200, 0, 1500,
                   80, 7, 1500, 0)*1
distances <- matrix(data = distances_vec, nr = 4, nc = 4)
avg_age <- c(40,20,50,25)
avg_house <- c(8,3,5,9)

# define parameters
beta_not <- 0.0004*(1/3.5)  # do we need to change this given we are moving from half weeks to days for each time step?  **************
beta_d <- 0.77*(1/3.5)  # do we need to change this given we are moving from half weeks to days for each time step?  **************
beta_ds <- 0
school <- 0
mu <- 0.23*(1/3.5)  # do we need to change this given we are moving from half weeks to days for each time step?  **************
nu <- 0
epsilon <- 1  # Stephen recommended sticking to 1 for now
n_locations <- 4
# infection_threshold <- 0.8   # not sure what this should be - ask Stephen for his thoughts
ro <- 96  # do we need to change this given we are moving from half weeks to days for each time step?  **************
beta_c <- 0.2  # based on data in papers linked here: https://docs.google.com/document/d/1MY5DfR6cU0gQ5wiKfxd1QaooSJZ4Io38A0uYwDPRO3U/edit
beta_a <- 0.2  # ******  maybe revisit this value ***********
gamma <- 0.125  # based on data in papers linked here: https://docs.google.com/document/d/1MY5DfR6cU0gQ5wiKfxd1QaooSJZ4Io38A0uYwDPRO3U/edit
                    # ~8 days until recovery
xi <- 0.1  # **** revisit this value ****
eta <- 0.1  # **** revisit this value ****

seed <- 0.1  # proportion infected to add at first infection step

# user inputs
init_Sc <- c(1, 1, .99, 1)
init_Ic <- c(0, 0, .01, 0)
init_Rc <- c(0, 0, 0, 0)

init_Sa <- c(1, 1, 1, 1)
init_Ia <- c(0, 0, 0, 0)
init_Ra <- c(0, 0, 0, 0)

t <- seq(from=0, to=100, by=1)  # change to from 100 to see if just taking a long time

model <- sir$new(beta_not=beta_not,
                 beta_d=beta_d,
                 beta_ds=beta_ds,
                 school=school,
                 mu=mu, nu=nu,
                 epsilon=epsilon,
                 n_locations=n_locations,
                 ro=ro,
                 beta_c=beta_c,
                 beta_a=beta_a,
                 gamma=gamma,
                 xi=xi,
                 eta=eta,
                 init_Sc=init_Sc,
                 init_Ic=init_Ic,
                 init_Rc=init_Rc,
                 init_Sa=init_Sa,
                 init_Ia=init_Ia,
                 init_Ra=init_Ra,
                 N=populations,
                 dij=distances,
                 avg_house=avg_house,
                 avg_age=avg_age,
                 seed=seed)

sol <- model$run(t)


# check if values add to 1
one_check <- sol %>%
  as.data.frame() %>%
  rename('Sc1' = 'Sc[1]',
         'Ic1' = 'Ic[1]',
         'Rc1' = 'Rc[1]',
         'Sa1' = 'Sa[1]',
         'Ia1' = 'Ia[1]',
         'Ra1' = 'Ra[1]') %>%
  mutate(c_one = Sc1 + Ic1 + Rc1, a_one = Sa1 + Ia1 + Ra1)

one_check$c_one
one_check$a_one  # ***** this doesn't add up to 1 yet *****


# plotting

sol_to_plot <- as_tibble(data.frame(model$run(t)))

# Generate a figure of the output: 
fig_sir <- sol_to_plot %>% 
  pivot_longer(names_to='population', values_to="n", cols=-c("step")) # %>% separate(col="population", into=c("Status","Population"))

ggplot(data=fig_sir, mapping=aes(x = step, y = n, col = population)) +
  geom_line() + theme(legend.position = "none")

fig_sir %>% filter(Population %in% c("1","2","3","4","5","6","7","8")) %>% 
  ggplot(aes(x = step, y = n, linetype=Status, col = Population)) +
  geom_line() + facet_wrap(~Population)






# Setting it up to run with Stephen's data


# read in data
coords <- as.matrix(read_csv("data/coords.csv", col_names = FALSE))
populations <- read_csv("data/pops09.csv", col_names = FALSE)$X1

# calculate great circle distances
distances <- as.matrix(gcd.slc(coords))

# set input parameters for Stephen's dataset
n_locations <- length(populations)

init_Sc <- rep(1,834)
init_Ic <- rep(0,834)
init_Rc <- rep(0,834)

sites_to_seed <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,100,200,300,400)

for (j in sites_to_seed){
  init_Sc[j] <- 0.99
  init_Ic[j] <- 0.01
}

# init_Sc[1] <- 0.99
# init_Ic[1] <- 0.01
