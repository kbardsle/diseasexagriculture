# --------------------------------------------------------
# calculating 3-digit zip population sizes and weighted latitude longitudes for 2017 census data
# 26 Apr 2024
# Katie Bardsley
# --------------------------------------------------------

library(tidyverse)
library(zipcodeR)

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






# --------------------------------------------------------
# calculating 3-digit zip population sizes and weighted latitude longitudes for 2017 census data
# 26 Apr 2024
# Katie Bardsley
# --------------------------------------------------------

library(tidyverse)
library(zipcodeR)  # to get state zip codes
  # some components being retired (or already have been?)
library(tigris)  # for state data
library(sf)

# get path to the script directory
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# Set the working directory to the parent directory of the script directory
setwd(file.path(script_dir, ".."))

# input state info
state_df <- read_csv("data/state_abbrev.csv") %>% filter(Contig == "Y")

zips_list <- lapply(state_df$Abbreviation, search_state)
zips_df <- bind_rows(zips_list)
zips_df <- zips_df %>% select(c("zipcode","state"))
colnames(zips_df)[1] <- "GEOID"

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

# add state info
zipdata_full_state <- full_join(zipdata_full, zips_df, by="GEOID")

# add column for 3 digit zip
zipdata_full_state <- zipdata_full_state %>%
  dplyr::mutate(ZIP3 = substr(as.character(GEOID), 1, 3))

# sum populations for 3-digit zips
zip3_pops <- zipdata_full_state %>% group_by(ZIP3) %>% summarize(POP3 = sum(POPULATION))

# combine with full df
zipdata_full3 <- full_join(zipdata_full_state, zip3_pops, by="ZIP3")

# filter out zips with less than 20,000 people
zipdata3_filtered <- zipdata_full3 %>% filter(POP3 >= 20000)

# calculate population-weighted lat/long
weighted_latlong <- zipdata3_filtered %>% 
  group_by(ZIP3) %>% 
  summarize(WEIGHTED_LAT = weighted.mean(x=INTPTLAT, w=POPULATION),
            WEIGHTED_LONG = weighted.mean(x=INTPTLONG, w=POPULATION))

# combine with full df
zipdata3_weighted_df <- full_join(zipdata3_filtered, weighted_latlong, by="ZIP3")
df_state_end <- zipdata3_weighted_df %>% 
  group_by(ZIP3) %>% 
  summarize(STATE = median(state), POP3 = mode(POP3), WEIGHTED_LAT = mode(WEIGHTED_LAT), WEIGHTED_LONG = mode(WEIGHTED_LONG))



df_state_end <- zipdata3_weighted_df %>%
  select(c("ZIP3","state","POP3","WEIGHTED_LAT", "WEIGHTED_LONG")) #%>% 
  #select()# %>% 
  distinct(ZIP3, state)


# simplified df with just 3-digit areas - doesn't have state info
# df_3digit <- full_join(zip3_pops, weighted_latlong, by="ZIP3")
# df_3digit_filtered <- df_3digit %>% filter(POP3 >= 20000)

# filter for contiguous US
df_US <- df_3digit_filtered %>% filter(WEIGHTED_LONG > -126 & WEIGHTED_LONG < -65) %>% filter(WEIGHTED_LAT > 23 & WEIGHTED_LAT < 50)

# save as csv
write.csv(df_US, "data/2017_pop_lat_long_data.csv", row.names=FALSE)
test <- read_csv("data/2017_pop_lat_long_data.csv")
