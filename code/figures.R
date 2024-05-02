# code for figures

# load packages
library(tidyverse)
library(wesanderson)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

# read in data (alternatively, get from disease model output directly)
infection_data <- read.csv("data/infection_start_data.csv")  # infection start timepoints
output_state_info_df <- read.csv("data/model_output_state_demographics.csv")  # also works where calls for fig_data or data_2017
state_data_raw <- read_csv("data/model_output_grouped_by_state.csv") # state level data
ag_demographic_data <- read_csv("data/migrants_merger.csv")  # demographic data for agricultural workforce
gen_demographic_data<- read_csv("data/general_population_demographics.csv")  # demographic data for general population


# MAP STATES TO REGIONS -------------------------------------

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


# PLOTTING -------------------------------------


# MAP OF SPATIAL INFECTION SPREAD -------------------------------------

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

ggsave("figures/infection_map.png", plot=map, width=8, height=5)


# DEMOGRAPHIC DATA HISTOGRAMS -------------------------------------

ag_clean_demo_data <- ag_demographic_data %>%
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

# normalize demographic data
gen_demographic_data_normalized <- gen_demographic_data %>%
  mutate(proportion_crowded = percent_crowded/100,
         proportion_w_kids = percent_with_children/100)

# histogram for crowding
hist_crowding <- ggplot() +
  geom_histogram(data = ag_clean_demo_data, aes(x = proportion_crowded, fill = "Agricultural Workforce"), alpha=0.6) +
  geom_histogram(data = gen_demographic_data_normalized, aes(x = proportion_crowded, fill = "General Population"), alpha=0.6) +
  scale_fill_manual(values = c("#4fc375", "#4574ca"), name = "Community") +
  labs(x = "proportion crowded") +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

ggsave("figures/prop_crowded_hist.png", plot=hist_crowding, width=5, height=3)

# histogram for households with children
hist_children <- ggplot() +
  geom_histogram(data = ag_clean_demo_data, aes(x = proportion_w_kids, fill = "Agricultural Workforce"), alpha=0.6) +
  geom_histogram(data = gen_demographic_data_normalized, aes(x = proportion_w_kids, fill = "General Population"), alpha=0.6) +
  scale_fill_manual(values = c("#4fc375", "#4574ca"), name = "Community") +
  labs(x = "proportion with children") +
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))

ggsave("figures/prop_children_hist.png", plot=hist_children, width=5, height=3)


# MAPS OF AGRICULTURAL WORKER DEMOGRAPHICS -------------------------------------

map_df <- output_state_info_df %>%  #data_2017 %>% 
  group_by(State) %>% 
  summarize(crowded=mean(proportion_crowded),
            children=mean(proportion_w_kids),
            state_ag_pop=mean(state_ag_population),
            ag_prop=(median(state_ag_population)/sum(POP3)))
map_df$State <- tolower(map_df$State)
colnames(map_df)[1] <- "region"

full_map_df <- full_join(us_map, map_df, by="region")

# map for household crowding
map_crowding <- ggplot() + 
  geom_polygon(data = full_map_df,
               aes(x = long, y = lat, group = group, fill = crowded)) +
  geom_polygon(data = map_data("state"), 
               aes(x = long, y = lat, group = group), 
               color = "black", fill = NA) +
  scale_fill_viridis() + 
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white"),
                          plot.background = element_rect(fill = "white"),
                          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  labs(fill="Proportion", title="Proportion of households with more than one resident per room")

ggsave("figures/crowding_map.png", plot=map_crowding, width=6, height=4)

# map for proportion of households with children
map_children <- ggplot() + 
  geom_polygon(data = full_map_df,
               aes(x = long, y = lat, group = group, fill = children)) +
  geom_polygon(data = map_data("state"), 
               aes(x = long, y = lat, group = group), 
               color = "black", fill = NA) +
  scale_fill_viridis() + 
  theme_minimal() + 
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  labs(fill="Proportion", title="Proportion of households with children")

ggsave("figures/children_map.png", plot=map_children, width=6, height=4)

# map for distribution of agricultural workers
map_ag_workers <- ggplot() + 
  geom_polygon(data = full_map_df,
               aes(x = long, y = lat, group = group, fill = ag_prop)) +
  geom_polygon(data = map_data("state"), 
               aes(x = long, y = lat, group = group), 
               color = "black", fill = NA) +
  scale_fill_viridis() + 
  theme_minimal() +
  labs(fill="Proportion", title="Proportion of state population that are agricultural workers")

ggsave("figures/ag_workers_map.png", plot=map_ag_workers, width=6, height=4)


# PEAK INFECTIOUS BY STATE, COLORED BY REGION -------------------------------------

# pull out peak infected value for general population
gen_peak_infected <- median(filter(state_data_clean, Community == "c")$peak_infected)

# make palette for plot
pal <- wes_palette("IsleofDogs1")
pal2 <- wes_palette("AsteroidCity2")
pal3 <- wes_palette("Darjeeling1")

peak_infect_bar <- state_data_clean %>%
  filter(Community == "a") %>%
  #arrange(peak_infected) %>%
  ggplot(aes(x = reorder(State_Abbreviation, -peak_infected), y = peak_infected, fill = Region)) +
  geom_col() +
  # scale_fill_viridis(discrete=TRUE, option="inferno") + 
  # scale_fill_manual(values = pal3) +
  scale_fill_brewer(palette = "Dark2") +
  geom_hline(yintercept = gen_peak_infected, linetype = 2, linewidth = 2, col = "black") +
  labs(x = "State", y = "Peak Proportion Infected") +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16))

ggsave("figures/peak_infection_state_bar.png", peak_infect_bar, width=14, height=7)


# SIR PLOT BY STATE (INFECTIOUS ONLY), COLORED BY REGION -------------------------------------

# plot infections curves

SIR_I_state <- state_data_raw %>% 
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
  filter(Disease_Status == "I") %>%  
  ggplot(aes(x=step, y=n_state, color=Region, linetype=Community, group=interaction(Region, State, Community))) + 
  geom_line() +
  # scale_color_viridis(discrete = TRUE, option = "turbo") +
  scale_color_brewer(palette = "Dark2") +
  scale_linetype_manual(labels = c("Agricultural Workers", "General Population"), values = c(1, 2)) +
  #scale_color_manual(values = reds_pal(50)) +
  theme_bw() + 
  labs(y="Proportion Infected", x="Timestep") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))

ggsave("figures/infectious_by_state.png", SIR_I_state, width=10, height=5)








# OLD FIGURES -------------------------------------

data_reformat <- output_state_info_df %>% pivot_wider(names_from = "Disease_Status", values_from = "n")
pops_4 <- sample(unique(data_reformat$Population), 4)
pops_50 <- sample(unique(data_reformat$Population), 50)

# figure with all pops plotted in the same figure, color by disease status
output_state_info_df %>% 
  ggplot(aes(x = step, y = n, col = Disease_Status, linetype=Community, group_by=Population)) +
  geom_line() + theme_bw()

# SIR plot 
# make palettes
reds_pal <- scales::hue_pal(h = c(0, 60))
greens_pal <- scales::hue_pal(h = c(80, 140))
blues_pal <- scales::hue_pal(h = c(200, 260))
scales::show_col(blues_pal(50))

# blue for s, red for I, green for R
SIR_no_legend <- data_reformat %>% filter(Population %in% pops_50) %>% ggplot() + 
  # geom_line(aes(x=step, y=S, color=Population, linetype=Community)) +
  #scale_color_stepsn(colors = c("#f5b0b0", "#400202")) +
  scale_color_manual(values = blues_pal(50)) +
  #scale_color_viridis(discrete = TRUE, option = "viridis") +
  new_scale_color() +
  geom_line(aes(x=step, y=I, color=Population, linetype=Community)) +
  # #scale_color_viridis(discrete = TRUE, option = "rocket") +
  # scale_color_manual(values = reds_pal(50)) +
  # new_scale_color() +
  # geom_line(aes(x=step, y=R, color=Population, linetype=Community)) + 
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


# figure with facets for populations
output_state_info_df %>% filter(Population %in% pops_4) %>% 
  ggplot(aes(x = step, y = n, col = Disease_Status, linetype=Community)) +
  geom_line() + facet_wrap(~Population) + theme_bw()

# SIR by state, colored by state
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