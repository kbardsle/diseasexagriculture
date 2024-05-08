# converting to hours lost from model output
# 5/8/24

library(tidyverse)

# get path to the script directory
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# Set the working directory to the parent directory of the script directory
setwd(file.path(script_dir, ".."))

# read in data
model_output_df <- read_csv("data/model_output_grouped_by_state.csv") %>%
  filter(Disease_Status == "I") # state level data

# add column for workers infected in a given timestep
model_output_df$workers_infected_timestep <- model_output_df$n_state*model_output_df$state_ag_pop


# map states to regions
region_mapping <- tibble(
  REGION6 = 1:6,
  RegionName = c("EAST", "SOUTHEAST", "MIDWEST", "SOUTHWEST", "NORTHWEST", "CALIFORNIA"),
  States = c(
    "North Carolina, Virginia, Kentucky, Tennessee, West Virginia, Connecticut, Maine, Massachusetts, New Hampshire, New York, Rhode Island, Vermont, Delaware, Maryland, New Jersey, Pennsylvania",
    "Arkansas, Louisiana, Mississippi, Alabama, Georgia, South Carolina, Florida",
    "Illinois, Indiana, Ohio, Iowa, Missouri, Kansas, Nebraska, North Dakota, South Dakota, Michigan, Minnesota, Wisconsin",
    "Arizona, New Mexico, Oklahoma, Texas",
    "Idaho, Montana, Wyoming, Colorado, Nevada, Utah, Oregon, Washington",
    "California"
  )
) %>%
  mutate(States = str_split(States, ", "))

# wrangle data
model_output_region_df <- model_output_df %>%
  filter(Disease_Status == "I") %>%
  mutate(State = str_to_title(State),
         Region = str_to_title(case_when(
           State %in% region_mapping[1,3][[1]][[1]] ~ region_mapping[1,2][[1]],
           State %in% region_mapping[2,3][[1]][[1]] ~ region_mapping[2,2][[1]],
           State %in% region_mapping[3,3][[1]][[1]] ~ region_mapping[3,2][[1]],
           State %in% region_mapping[4,3][[1]][[1]] ~ region_mapping[4,2][[1]],
           State %in% region_mapping[5,3][[1]][[1]] ~ region_mapping[5,2][[1]],
           State %in% region_mapping[6,3][[1]][[1]] ~ region_mapping[6,2][[1]],
           TRUE ~ NA
         )))

# add column for average hours worked per week (based on regional info)
model_output_region_df <- model_output_region_df %>%
  mutate(avg_week_hours <- ifelse(Region == "East", 44.9,
                                 ifelse(Region == "Southeast", 42.4,
                                        ifelse(Region == "Midwest", 41.1,
                                               ifelse(Region == "Southwest", 48.3,
                                                      ifelse(Region == "Northwest", 42.8,
                                                             ifelse(Region == "California", 44.6, 0)))))))
model_output_region_df$avg_day_hours <- model_output_region_df$avg_week_hours/7

# calculate the hours lost per timestep
model_output_region_df$hours_lost_timestep <- model_output_region_df$workers_infected_timestep * model_output_region_df$avg_day_hours

# filter for just agricultural workforce
model_output_region_ag <- model_output_region_df %>% filter(Community == "a")

# add column for total hours
model_output_region_ag$total_state_hours <- model_output_region_ag$state_ag_pop * model_output_region_ag$avg_week_hours * 52

# sum up by state
state_hrs_lost <- model_output_region_ag %>% group_by(State) %>% summarize(hours_lost = sum(hours_lost_timestep),
                                                         hours_total = mean(total_state_hours), 
                                                         prop_yearly_hrs_lost = hours_lost/hours_total)

write.csv(state_hrs_lost, "data/state_hours_lost.csv", row.names=FALSE)
