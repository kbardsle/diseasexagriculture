# Defining spatially and temporally explicit disease model using Odin
# Group rotation project - Spring 2024

# LOAD LIBRARIES -------------------------------------

library(odin)
library(tidyverse)
library(codep)  # for calculating great circle distances
library(maps)  # for plotting the map
library(mapdata)  # for plotting the map
library(viridis)
library(ggnewscale)
library(cowplot)



# SET WORKING DIRECTORY -------------------------------------

# get path to the script directory
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# Set the working directory to the parent directory of the script directory
setwd(file.path(script_dir, ".."))


# BUILD MODEL -------------------------------------

# build spatial model
sir <- odin::odin({

  # Community transmission:
  
  # Logic for community-level SIR Calculations:
    # if a location is newly infected (infected[i] > 0 and I[i] == 0), 
        # increase the proportion of infected by a set seed (here 0.01) to seed the infection
    # if a location was infected in an earlier time step, calculate derivatives as normal
    # if a location is still not infected, set value equal to previous time step

  # SUSCEPTIBLE:
  update(Sc[]) <- if (infected[i] > 0 && Ic[i] == 0) Sc[i] - (beta_base+xi*proportion_crowded_c[i]+eta*proportion_w_kids_c[i])*Sc[i]*(Ic[i]+seed) else if (Ic[i] > 0) Sc[i] - (beta_base+xi*proportion_crowded_c[i]+eta*proportion_w_kids_c[i])*Sc[i]*Ic[i] else Sc[i]

  # INFECTED:
  update(Ic[]) <- if (infected[i] > 0 && Ic[i] == 0) Ic[i] + (beta_base+xi*proportion_crowded_c[i]+eta*proportion_w_kids_c[i])*Sc[i]*(Ic[i]+seed) - gamma*(Ic[i]+seed) else if (Ic[i] > 0 ) Ic[i] + (beta_base+xi*proportion_crowded_c[i]+eta*proportion_w_kids_c[i])*Sc[i]*Ic[i] - gamma*Ic[i] else Ic[i]
  
  # RECOVERED:
  update(Rc[]) <- if (infected[i] > 0 && Ic[i] == 0) Rc[i] + gamma*(Ic[i]+seed) else if (Ic[i] > 0 ) Rc[i] + gamma*Ic[i] else Rc[i]

  # Agricultural workforce transmission:
  
  # Logic for agricultural workforce SIR Calculations:
    # if the new calculated value is less than zero, set it to 0
    # if the new calculated value is greater than one, set it to 1
    # if the new calculated value is between 0 and 1 inclusive, 
        # use the calculated value as the new value
  
  # SUSCEPTIBLE:
  update(Sa[]) <- if (Sa[i] - (((beta_base+xi*proportion_crowded_c[i]+eta*proportion_w_kids_c[i])*resistance*Ic[i]*(1-assortment_prob))+((beta_base+xi*proportion_crowded_a[i]+eta*proportion_w_kids_a[i])*resistance*Ia[i]*assortment_prob))*Sa[i] < 0) 0 else if (Sa[i] - (((beta_base+xi*proportion_crowded_c[i]+eta*proportion_w_kids_c[i])*resistance*Ic[i]*(1-assortment_prob))+((beta_base+xi*proportion_crowded_a[i]+eta*proportion_w_kids_a[i])*resistance*Ia[i]*assortment_prob))*Sa[i] > 1) 1 else Sa[i] - (((beta_base+xi*proportion_crowded_c[i]+eta*proportion_w_kids_c[i])*resistance*Ic[i]*(1-assortment_prob))+((beta_base+xi*proportion_crowded_a[i]+eta*proportion_w_kids_a[i])*resistance*Ia[i]*assortment_prob))*Sa[i]
  
  # INFECTED:
  update(Ia[]) <- if (Ia[i] + (((beta_base+xi*proportion_crowded_c[i]+eta*proportion_w_kids_c[i])*resistance*Ic[i]*(1-assortment_prob))+((beta_base+xi*proportion_crowded_a[i]+eta*proportion_w_kids_a[i])*resistance*Ia[i]*assortment_prob))*Sa[i] - gamma*Ia[i] < 0) 0 else if (Ia[i] + (((beta_base+xi*proportion_crowded_c[i]+eta*proportion_w_kids_c[i])*resistance*Ic[i]*(1-assortment_prob))+((beta_base+xi*proportion_crowded_a[i]+eta*proportion_w_kids_a[i])*resistance*Ia[i]*assortment_prob))*Sa[i] - gamma*Ia[i] > 1) 1 else Ia[i] + (((beta_base+xi*proportion_crowded_c[i]+eta*proportion_w_kids_c[i])*resistance*Ic[i]*(1-assortment_prob))+((beta_base+xi*proportion_crowded_a[i]+eta*proportion_w_kids_a[i])*resistance*Ia[i]*assortment_prob))*Sa[i] - gamma*Ia[i]
  
  # RECOVERED:
  update(Ra[]) <- if (Ra[i] + gamma*Ia[i] < 0) 0 else if (Ra[i] + gamma*Ia[i] > 1) 1 else Ra[i] + gamma*Ia[i]
  
  # DEFINE TERMS:
  
  k_dij[,] <- exp(-dij[i,j]/ro)  # exponential kernel
  
  # matrix of the values for the equation numerator for all combinations of i and j:
  NK_num[,] <- N[j]^nu * k_dij[i,j] * (if (Ic[j] > 0) 1 else 0)
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
  lambda[] <- if (Ic[i] == 0) (beta_not + (beta_d + (school * beta_ds)) * N[i]^mu * (sum(NK_num[i,]))/(denom[i])) else 0
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
  # print("infection probs: {infection_prob[1]}, {infection_prob[2]}, {infection_prob[3]}, {infection_prob[4]}, {infection_prob[5]}, {infection_prob[6]}, {infection_prob[7]}, {infection_prob[805]}")
  # print("lambda: {lambda[1]}, {lambda[2]}, {lambda[3]}, {lambda[4]}")
  # print("dij: {dij[1,1]}, {dij[1,2]}, {dij[1,3]}, {dij[1,4]}")
  # print("k_dij: {k_dij[1,1]}, {k_dij[1,2]}, {k_dij[1,3]}, {k_dij[1,4]}")
  # print("NK_num: {NK_num[1,1]}, {NK_num[2,1]}, {NK_num[2,3]}, {NK_num[3,1]}, {NK_num[4,1]}, {NK_num[5,1]}")
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
  # beta_c <- user()
  beta_base <- user()
  gamma <- user()
  xi <- user()
  eta <- user()
  seed <- user()
  assortment_prob <- user(0.5)
  resistance <- user(1)
  
  # agricultural workforce community parameters
  proportion_w_kids_a[] <- user()
  proportion_crowded_a[] <- user()
  
  # general community parameters
  proportion_w_kids_c[] <- user()
  proportion_crowded_c[] <- user()
  
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
  dim(proportion_crowded_a) <- n_locations
  dim(proportion_w_kids_a) <- n_locations
  dim(proportion_crowded_c) <- n_locations
  dim(proportion_w_kids_c) <- n_locations
  
  # dim(x) <- n_locations  # added to help with debugging
  
}, debug_enable=TRUE)


# SET INPUT PARAMETERS & DATA -------------------------------------

# DEFINE PARAMETERS

# parameters for spatial model
beta_not <- 0.0004*(1/3.5)  # adjusted for half weeks --> days
beta_d <- 0.77*(1/3.5)  # adjusted for half weeks --> days
beta_ds <- 0
school <- 0
mu <- 0.23  # likely don't have to adjust for days/half weeks
nu <- 0
epsilon <- 1  # Stephen recommended sticking to 1 for now
ro <- 96


# beta_c <- 0.2  # based on data in papers linked here: 
                    # https://docs.google.com/document/d/1MY5DfR6cU0gQ5wiKfxd1QaooSJZ4Io38A0uYwDPRO3U/edit
beta_base <- 0.1 
gamma <- 0.125  # based on data in papers linked here: 
                    # https://docs.google.com/document/d/1MY5DfR6cU0gQ5wiKfxd1QaooSJZ4Io38A0uYwDPRO3U/edit
                    # ~8 days until recovery
xi <- 0.92  # coefficient for proportion crowded - based on parameterization
eta <- 0.25  # coefficient for proportion with kids - based on parameterization

seed <- 0.01  # proportion infected to add at first infection step

t <- seq(from=0, to=350, by=1)


# INPUT DATA

# use 2017 population and latitude/longitude data
data_2017 <- read_csv("data/2017_pop_demo_data_both_communities.csv")
coords <- as.matrix(data_2017 %>% select(c("WEIGHTED_LONG","WEIGHTED_LAT")))
populations <- data_2017$POP3

avg_pop <- mean(populations)
normalized_populations <- populations/avg_pop

# calculate great circle distances
distances <- as.matrix(gcd.slc(coords))
# distances <- distances * 1000

# set input parameters for 2017 data
n_locations <- length(populations)

init_Sc <- rep(1,n_locations)
init_Ic <- rep(0,n_locations)
init_Rc <- rep(0,n_locations)

init_Sa <- rep(1,n_locations)
init_Ia <- rep(0,n_locations)
init_Ra <- rep(0,n_locations)

sites_to_seed <- c(350, 359, 614, 835)
# Grenada, MS - 350
# Albany, GA - 359
# Stockton, CA - 614
# Omaha, NE - 835

for (j in sites_to_seed){
  init_Sc[j] <- 0.99
  init_Ic[j] <- 0.01
}


# RUN MODEL -------------------------------------

model <- sir$new(beta_not=beta_not,
                 beta_d=beta_d,
                 beta_ds=beta_ds,
                 school=school,
                 mu=mu, nu=nu,
                 epsilon=epsilon,
                 n_locations=n_locations,
                 ro=ro,
                 # beta_c=beta_c,
                 beta_base=beta_base,
                 gamma=gamma,
                 xi=xi,
                 eta=eta,
                 init_Sc=init_Sc,
                 init_Ic=init_Ic,
                 init_Rc=init_Rc,
                 init_Sa=init_Sa,
                 init_Ia=init_Ia,
                 init_Ra=init_Ra,
                 N=normalized_populations,
                 dij=distances,
                 proportion_crowded_a=data_2017$proportion_crowded_a,
                 proportion_w_kids_a=data_2017$proportion_w_kids_a,
                 proportion_crowded_c=data_2017$proportion_crowded_c,
                 proportion_w_kids_c=data_2017$proportion_w_kids_c,
                 seed=seed)

# sol <- model$run(t)

# FORMAT AND SAVE MODEL OUTPUT -------------------------------------

sol_to_plot <- as_tibble(data.frame(model$run(t)))

# make reformatted dataframe with separate columns for:
    # "Disease_Status" (S, I, or R)
    # "Community" (c for general community, a for agricultural workforce)
    # "Population" for population number
    # "step" for timestep
    # "n" for value
fig_data <- sol_to_plot %>% 
  pivot_longer(names_to='population', values_to="n", cols=-c("step"))  %>% separate(col="population", into=c("Status","Population")) %>% separate(col="Status", into=c("Disease_Status","Community"),sep=1:1)

# dataframe for infection start date
ind <- do.call(rbind, lapply(870:1737, function(i){
  data.frame(population = i-869, infection_start_index = min(which(sol_to_plot[,i] != 0)))
}))
infection_data <- data.frame(ind$population, coords, ind$infection_start_index) %>% 
  rename(latitude = WEIGHTED_LAT,
         longitude = WEIGHTED_LONG,
         infection_start_index = ind.infection_start_index,
         population = ind.population)

# save
out_name_infection_data <- paste0("data/infection_start_data_", format(Sys.Date(), "%m.%d.%y"), ".csv")
write.csv(infection_data, out_name_infection_data)

# add combined dataframe with state information
data_2017_temp <- data_2017
colnames(data_2017_temp)[13] <- "Population"
data_2017_temp$Population <- as.character(data_2017_temp$Population)

output_state_info_df <- full_join(fig_data, data_2017_temp, by="Population")

# save
out_name_state_info <- paste0("data/model_output_state_demographics_", format(Sys.Date(), "%m.%d.%y"), ".csv")
write.csv(output_state_info_df, out_name_state_info)

# state level data
state_output_df <- output_state_info_df %>% group_by(step, Disease_Status, Community, State, State_Abbreviation) %>% 
  summarize(state_ag_pop=median(state_ag_population),
            state_pop=sum(POP3),
            n_state=weighted.mean(n, POP3))

# save
out_name_state_output <- paste0("data/model_output_grouped_by_state_", format(Sys.Date(), "%m.%d.%y"), ".csv")
write.csv(state_output_df, out_name_state_output)



# TESTS -------------------------------------

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
one_check$a_one

# looking at SIR values
fig_data %>% group_by(Population, Disease_Status) %>% summarize(min=min(n), max=max(n), mean=mean(n))

# checking first time step for location infection
summary(ind)

ggplot(data=ind, mapping=aes(x=infection_start_index)) + geom_histogram()

# distance from springfield
spring_dist <- distances[,1]
data.frame(ind=ind$infection_start_index, dist=spring_dist) %>% 
  ggplot(mapping=aes(x=dist, y=ind)) + geom_point()

# population size 
data.frame(ind=ind$infection_start_index, pop=populations) %>% 
  ggplot(mapping=aes(x=pop, y=ind)) + geom_point()