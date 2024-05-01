# load libraries
library(tidyverse)
library(wesanderson)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

# read in state level data
state_data_raw <- read_csv("data/model_output_grouped_by_state.csv")

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
state_data_clean <- state_data_raw %>%
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
         ))) %>%
  group_by(Region, State, State_Abbreviation, Disease_Status, Community) %>%
  summarize(peak_infected = max(n_state))


# pull out peak infected value for general population
gen_peak_infected <- median(filter(state_data_clean, Community == "c")$peak_infected)

# make palette for plot
pal <- wes_palette("IsleofDogs1")
pal2 <- wes_palette("AsteroidCity2")

state_data_clean %>%
  filter(Community == "a") %>%
  #arrange(peak_infected) %>%
  ggplot(aes(x = reorder(State_Abbreviation, -peak_infected), y = peak_infected, fill = Region)) +
  geom_col() +
  scale_fill_manual(values = pal) +
  geom_hline(yintercept = gen_peak_infected, linetype = 2, linewidth = 3, col = "red") +
  labs(x = "State", y = "Peak Proportion Infected") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))

# plot infections curves
ggplot(data = filter(state_data_raw, Disease_Status == "I"), aes(x=step, y=n_state, color=State_Abbreviation, linetype=Community)) + 
  geom_line() +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  scale_linetype_manual(labels = c("Agricultural Workers", "General Population"), values = c(1, 2)) +
  #scale_color_manual(values = reds_pal(50)) +
  theme_bw() + 
  labs(y="Proportion Infected", x="Timestep") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))
