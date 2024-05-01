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
  update(Sc[]) <- if (infected[i] > 0 && Ic[i] == 0) Sc[i] - beta_c*Sc[i]*(Ic[i]+seed) else if (Ic[i] > 0) Sc[i] - beta_c*Sc[i]*Ic[i] else Sc[i]

  # INFECTED:
  update(Ic[]) <- if (infected[i] > 0 && Ic[i] == 0) Ic[i] + beta_c*Sc[i]*(Ic[i]+seed) - gamma*(Ic[i]+seed) else if (Ic[i] > 0 ) Ic[i] + beta_c*Sc[i]*Ic[i] - gamma*Ic[i] else Ic[i]
  
  # RECOVERED:
  update(Rc[]) <- if (infected[i] > 0 && Ic[i] == 0) Rc[i] + gamma*(Ic[i]+seed) else if (Ic[i] > 0 ) Rc[i] + gamma*Ic[i] else Rc[i]

  # Agricultural workforce transmission:
  
  # Logic for agricultural workforce SIR Calculations:
    # if the new calculated value is less than zero, set it to 0
    # if the new calculated value is greater than one, set it to 1
    # if the new calculated value is between 0 and 1 inclusive, 
        # use the calculated value as the new value
  
  # SUSCEPTIBLE:
  update(Sa[]) <- if (Sa[i] - ((beta_c*Ic[i]*(1-assortment_prob))+((beta_a+xi*proportion_crowded[i]+eta*proportion_w_kids[i])*Ia[i]*assortment_prob))*Sa[i] < 0) 0 else if (Sa[i] - ((beta_c*Ic[i]*(1-assortment_prob))+((beta_a+xi*proportion_crowded[i]+eta*proportion_w_kids[i])*Ia[i]*assortment_prob))*Sa[i] > 1) 1 else Sa[i] - ((beta_c*Ic[i]*(1-assortment_prob))+((beta_a+xi*proportion_crowded[i]+eta*proportion_w_kids[i])*Ia[i]*assortment_prob))*Sa[i]
  
  # INFECTED:
  update(Ia[]) <- if (Ia[i] + ((beta_c*Ic[i]*(1-assortment_prob))+((beta_a+xi*proportion_crowded[i]+eta*(1-proportion_w_kids[i]))*Ia[i]*assortment_prob))*Sa[i] - gamma*Ia[i] < 0) 0 else if (Ia[i] + ((beta_c*Ic[i]*(1-assortment_prob))+((beta_a+xi*proportion_crowded[i]+eta*(1-proportion_w_kids[i]))*Ia[i]*assortment_prob))*Sa[i] - gamma*Ia[i] > 1) 1 else Ia[i] + ((beta_c*Ic[i]*(1-assortment_prob))+((beta_a+xi*proportion_crowded[i]+eta*(1-proportion_w_kids[i]))*Ia[i]*assortment_prob))*Sa[i] - gamma*Ia[i]
  
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
  beta_c <- user()
  beta_a <- user()
  gamma <- user()
  xi <- user()
  eta <- user()
  seed <- user()
  assortment_prob <- user(0.5)
  
  # agricultural workforce community parameters
  proportion_w_kids[] <- user()
  proportion_crowded[] <- user()
  
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
  dim(proportion_crowded) <- n_locations
  dim(proportion_w_kids) <- n_locations
  
  # dim(x) <- n_locations  # added to help with debugging
  
}, debug_enable=TRUE)


# SET INPUT PARAMETERS & DATA -------------------------------------

# DEFINE PARAMETERS

beta_not <- 0.0004*(1/3.5)  # adjusted for half weeks --> days
beta_d <- 0.77*(1/3.5)  # adjusted for half weeks --> days
beta_ds <- 0
school <- 0
mu <- 0.23  # likely don't have to adjust for days/half weeks
nu <- 0
epsilon <- 1  # Stephen recommended sticking to 1 for now
ro <- 96
beta_c <- 0.2  # based on data in papers linked here: 
                    # https://docs.google.com/document/d/1MY5DfR6cU0gQ5wiKfxd1QaooSJZ4Io38A0uYwDPRO3U/edit
beta_a <- 0.1  # ******  maybe revisit this value ***********
gamma <- 0.125  # based on data in papers linked here: 
                    # https://docs.google.com/document/d/1MY5DfR6cU0gQ5wiKfxd1QaooSJZ4Io38A0uYwDPRO3U/edit
                    # ~8 days until recovery
xi <- 0.92  # coefficient for proportion crowded **** revisit this value after parameter fitting ****
eta <- 0.25  # coefficient for proportion with kids **** revisit this value after parameter fitting ****

seed <- 0.01  # proportion infected to add at first infection step

t <- seq(from=0, to=150, by=1)


# INPUT DATA

# use basic test data set - 4 locations

populations <- c(100, 100, 100, 100)
normalized_populations <- c(1,1,1,1)
distances_vec <- c(0, 5, 10, 80,
                   5, 0, 200, 7,
                   10, 200, 0, 1500,
                   80, 7, 1500, 0)*1
distances <- matrix(data = distances_vec, nr = 4, nc = 4)
proportion_w_kids <- c(0.3,0.5,0.7,0.9)
proportion_crowded <- c(0.3,0.5,0.7,0.9)
n_locations <- 4

init_Sc <- c(1, 1, .99, 1)
init_Ic <- c(0, 0, .01, 0)
init_Rc <- c(0, 0, 0, 0)

init_Sa <- c(1, 1, 1, 1)
init_Ia <- c(0, 0, 0, 0)
init_Ra <- c(0, 0, 0, 0)

# use test data from Stephen

# read in data
coords <- as.matrix(read_csv("data/coords.csv", col_names = FALSE))
populations <- read_csv("data/pops09.csv", col_names = FALSE)$X1
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




avg_pop <- mean(populations)
normalized_populations <- populations/avg_pop

# calculate great circle distances
distances <- as.matrix(gcd.slc(coords))
# distances <- distances * 1000

# set input parameters for Stephen's dataset
n_locations <- length(populations)

init_Sc <- rep(1,834)
init_Ic <- rep(0,834)
init_Rc <- rep(0,834)

init_Sa <- rep(1,834)
init_Ia <- rep(0,834)
init_Ra <- rep(0,834)

sites_to_seed <- c(1)  #,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,100,200,300,400)

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
                 N=normalized_populations,
                 dij=distances,
                 proportion_crowded=clean_demo_data$proportion_crowded,
                 proportion_w_kids=clean_demo_data$proportion_w_kids,
                 seed=seed)

sol <- model$run(t)

# FORMAT MODEL OUTPUT -------------------------------------

sol_to_plot <- as_tibble(data.frame(model$run(t)))

# make reformatted dataframe with separate columns for:
    # "Disease_Status" (S, I, or R)
    # "Community" (c for general community, a for agricultural workforce)
    # "Population" for population number
    # "step" for timestep
    # "n" for value
fig_data <- sol_to_plot %>% 
  pivot_longer(names_to='population', values_to="n", cols=-c("step"))  %>% separate(col="population", into=c("Status","Population")) %>% separate(col="Status", into=c("Disease_Status","Community"),sep=1:1)

# make reformatted dataframe with separate columns for S, I, and R:
data_reformat <- fig_data %>% pivot_wider(names_from = "Disease_Status", values_from = "n")


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
ind <- do.call(rbind, lapply(836:1669, function(i){
  data.frame(population = i-835, infection_start_index = min(which(sol_to_plot[,i] != 0)))
}))

summary(ind)

ggplot(data=ind, mapping=aes(x=infection_start_index)) + geom_histogram()

# distance from springfield
spring_dist <- distances[,1]
data.frame(ind=ind$infection_start_index, dist=spring_dist) %>% 
  ggplot(mapping=aes(x=dist, y=ind)) + geom_point()

# population size 
data.frame(ind=ind$infection_start_index, pop=populations) %>% 
  ggplot(mapping=aes(x=pop, y=ind)) + geom_point()


# PLOTTING -------------------------------------

pops_50 <- sample(unique(data_reformat$Population), 50)
pops_4 <- sample(unique(data_reformat$Population), 4)

# figure with subsample of 50 pops with different colors for each population
data_reformat %>% filter(Population %in% pops_50) %>% 
ggplot() + 
  geom_line(aes(x=step, y=S, color=Population, linetype=Community)) + 
  geom_line(aes(x=step, y=I, color=Population, linetype=Community)) + 
  geom_line(aes(x=step, y=R, color=Population, linetype=Community)) + 
  theme_bw() + guides(color=FALSE)

# ggplot(data=fig_sir, mapping=aes(x = step, y = n, col = population)) +
  # geom_line() #+ theme(legend.position = "none")

# figure with facets for populations
fig_data %>% filter(Population %in% pops_4) %>% 
  ggplot(aes(x = step, y = n, col = Disease_Status, linetype=Community)) +
  geom_line() + facet_wrap(~Population) + theme_bw()

# figure with all pops plotted in the same figure, color by disease status
fig_data %>% 
  ggplot(aes(x = step, y = n, col = Disease_Status, linetype=Community, group_by=Population)) +
  geom_line() + theme_bw()

# map plot
infection_data <- data.frame(coords, ind$infection_start_index)

infection_data <- infection_data %>% rename(latitude = X2,
       longitude = X1,
       infection_start_index = ind.infection_start_index)

us_map <- map_data("state")

map <- ggplot() +
  geom_polygon(data = us_map, 
               aes(x = long, y = lat, group = group), 
               fill = "white") +
  geom_polygon(data = map_data("state"), 
               aes(x = long, y = lat, group = group), 
               color = "black", fill = NA) +
  geom_point(data = infection_data, 
             aes(x = longitude, y = latitude, color = infection_start_index), 
             size = 1.5, alpha=0.85) + 
  scale_color_viridis(direction = -1, name = "Days to Infection") + 
  labs(title = "Infection spread across locations") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),  # Set plot background to white
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))  # Set plot margins to ensure white background around the plot
map
# ggsave("infection_map_4.22.24.png", plot=map, width=8, height=5)


# SIR plot 
# make palettes
reds_pal <- scales::hue_pal(h = c(0, 60))
greens_pal <- scales::hue_pal(h = c(80, 140))
blues_pal <- scales::hue_pal(h = c(200, 260))
scales::show_col(blues_pal(50))

# blue for s, red for I, green for R
SIR_no_legend <- data_reformat %>% filter(Population %in% pops_50) %>% ggplot() + 
  geom_line(aes(x=step, y=S, color=Population, linetype=Community)) +
  #scale_color_stepsn(colors = c("#f5b0b0", "#400202")) +
  scale_color_manual(values = blues_pal(50)) +
  #scale_color_viridis(discrete = TRUE, option = "viridis") +
  new_scale_color() +
  geom_line(aes(x=step, y=I, color=Population, linetype=Community)) +
  #scale_color_viridis(discrete = TRUE, option = "rocket") +
  scale_color_manual(values = reds_pal(50)) +
  new_scale_color() +
  geom_line(aes(x=step, y=R, color=Population, linetype=Community)) + 
  #scale_color_viridis(discrete = TRUE, option = "mako") +
  scale_color_manual(values = greens_pal(50)) +
  theme_bw() + 
  labs(y="Population proportion", x="Timestep") +
  guides(color=FALSE) +
  theme(legend.position = "none")

# make legend
legend_plot <- get_legend(fig_data %>% 
  filter(Population %in% pops_50) %>% 
  mutate(Disease_Status = fct_relevel(Disease_Status, c("S", "I", "R"))) %>%
  ggplot(aes(x = step, y = n, color = Disease_Status, linetype = Community)) + 
  geom_line() +
  scale_color_manual(values = c(S = "blue", I = "red", R = "green")) +
  scale_linetype_manual(labels = c("Agricultural Workers", "General Population"), values = c(1, 2)) +
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 15))
  )

plot_grid(SIR_no_legend, legend_plot, rel_widths = c(0.8, 0.2))


# PREPARE TO INPUT AG DATA -------------------------------------

###################################
# FUNCTION: normalize_ag_data
# description: normalize data in column to be between 0 and 1
# inputs: column to be normalized
# outputs: normalized column
# ---------------------------------
normalize_ag_data <- function(input_col){
  
  min <- min(input_col)
  range <- max(input_col) - min
  
  output_col <- (input_col - min)/range
  
  return(output_col)

}
###################################







