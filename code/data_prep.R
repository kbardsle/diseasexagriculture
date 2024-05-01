# data prep for use in disease model

# load packages
library(tidyverse)
library(zipcodeR)

# get path to the script directory
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# Set the working directory to the parent directory of the script directory
setwd(file.path(script_dir, ".."))

# calculating 3-digit zip population sizes for 2017 census data -------------------------------------

# read in data and reformat
latlong_zip <-read.table('data/2017_Gaz_zcta_national.txt',
                         sep='\t',
                         header=TRUE,
                         colClasses = "character")
# data from: 
# https://www.census.gov/geographies/reference-files/time-series/geo/gazetteer-files.2017.html#list-tab-264479560
# see here for more information on column contents (very bottom of page): 
# https://www.census.gov/programs-surveys/geography/technical-documentation/records-layout/gaz-record-layouts.2017.html#list-tab-1913338080

pop_zip <- read_csv("data/ACSDP5Y2017.DP05-Data.csv")

pop_zip_simple <- pop_zip %>% slice(-c(1)) %>% separate(col="NAME", into=c("prefix","GEOID")) %>% select(c(GEOID, DP05_0001E))
colnames(pop_zip_simple)[2] <- "POPULATION"

# combine dataframes
zipdata_full <- full_join(pop_zip_simple, latlong_zip, by="GEOID")

# convert lat long, pop to numeric
zipdata_full$INTPTLAT <- as.numeric(zipdata_full$INTPTLAT)
zipdata_full$INTPTLONG <- as.numeric(zipdata_full$INTPTLONG)
zipdata_full$POPULATION <- as.numeric(zipdata_full$POPULATION)

# add column for 3 digit zip
zipdata_full <- zipdata_full %>%
  dplyr::mutate(ZIP3 = substr(as.character(GEOID), 1, 3))

# sum populations for 3-digit zips
zip3_pops <- zipdata_full %>% group_by(ZIP3) %>% summarize(POP3 = sum(POPULATION))

# combine with full df
zipdata_full3 <- full_join(zipdata_full, zip3_pops, by="ZIP3")

# filter out zips with less than 20,000 people
zipdata3_filtered <- zipdata_full3 %>% filter(POP3 >= 20000)

# calculating 3-digit zip code population weighted lat long -------------------------------------

# calculate population-weighted lat/long
weighted_latlong <- zipdata3_filtered %>% 
  group_by(ZIP3) %>% 
  summarize(WEIGHTED_LAT = weighted.mean(x=INTPTLAT, w=POPULATION),
            WEIGHTED_LONG = weighted.mean(x=INTPTLONG, w=POPULATION))

# combine with full df
zipdata3_weighted_df <- full_join(zipdata3_filtered, weighted_latlong, by="ZIP3")

# simplified df with just 3-digit areas
df_3digit <- full_join(zip3_pops, weighted_latlong, by="ZIP3")
df_3digit_filtered <- df_3digit %>% filter(POP3 >= 20000)

# filter for contiguous US
df_US <- df_3digit_filtered %>% filter(WEIGHTED_LONG > -126 & WEIGHTED_LONG < -65) %>% filter(WEIGHTED_LAT > 23 & WEIGHTED_LAT < 50)

# adding state information to 3-digit zip dataset -------------------------------------

# input state info
state_df <- read_csv("data/state_abbrev.csv") %>% filter(Contig == "Y")

zips_list <- lapply(state_df$Abbreviation, search_state)
zips_df <- bind_rows(zips_list)
zips_df <- zips_df %>% select(c("zipcode","state"))

zips_df <- zips_df %>%
  dplyr::mutate(ZIP3 = substr(as.character(zipcode), 1, 3)) %>% 
  select(c("state","ZIP3")) %>% distinct()

# bind state info with pop/lat long info
df_US_states <- full_join(df_US, zips_df, by="ZIP3") %>% filter(!is.na(POP3))

# save as csv
write.csv(df_US_states, "data/2017_pop_lat_long_data_states.csv", row.names=FALSE)
# df_US_states <- read.csv("data/2017_pop_lat_long_data_states.csv")

# combine with agricultural worker demographic data -------------------------------------

# demographic data from: https://github.com/naiacasina/migrants_R_proj

# combine with demographic data for ag workers
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
            proportion_w_kids = weighted.mean(proportion_w_kids, Value),
            state_ag_population = sum(Value))

# add state names
colnames(state_df)[2] <- "State_Abbreviation"
colnames(df_US_states)[5] <- "State_Abbreviation"
df_US_states_names <- full_join(df_US_states, state_df, by="State_Abbreviation")
df_US_states_names$State <- toupper(df_US_states_names$State)

# do join with demographic data
data_complete_ag <- full_join(df_US_states_names, clean_demo_data, by="State")

# filter out washington dc
data_complete_ag <- data_complete_ag %>% filter(State_Abbreviation != "DC")

# add population ID
data_complete_ag <- data_complete_ag %>% 
  mutate(ID = seq.int(nrow(data_complete_ag)))

# save as csv - input for disease model
write.csv(data_complete_ag, "data/2017_pop_demo_data_agricultural.csv", row.names=FALSE)




# figuring out which locations to seed

# identifying location indexes for 4 sites from Stephen's paper
# Grenada, MS - 389
# Albany, GA - 398
# Stockton, CA - 952
# Omaha, NE - 681

seed_zips <- c(389, 398, 952, 681)

seed_site_df <- data_complete_ag %>% 
  filter(ZIP3 %in% seed_zips)

seed_indices <- seed_site_df$ID
