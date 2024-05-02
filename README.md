# diseasexagriculture
Repository for IQ Bio group rotation.

## Project Abstract

**Assessing the impact of an emerging respiratory pathogen on US agricultural workers and food production** 

Agriculture is a critical component of the US economy, contributing over 200 billion dollars annually to GDP. Many crops rely heavily on manual labor, but the impact of disruptions to the agricultural workforce on agricultural productivity is not well understood. Here we study the effects of the outbreak of a novel respiratory pathogen on the US agricultural workforce and its cascading impacts on agricultural output and economic production. We investigate how demographic variables, including household crowding and age structure, of the agricultural workforce impact disease transmission by combining a spatially and temporally explicit disease model with two separate SIR (susceptible-infectious-recovered) compartmental models of disease dynamics. This multi-layer approach makes it possible to model the patterns of disease transmission across the United States as well as within-community transmission for both the general population and the agricultural workforce at the state scale. Our model is parameterized using data from the 2017 US Census and National Agricultural Workers Survey (NAWS) from the same year. This model can be used to explore the effects of a novel disease outbreak on the US agricultural workforce and the economic impacts of the resulting agricultural productivity losses. By studying disruptions to this workforce, policies can be designed for both cost-effective and targeted disease surveillance strategies, as well as public health interventions specifically tailored to vulnerabilities and disease spread for this population. Additionally, this model will guide economic policy surrounding expected changes in food supply and the impacts on GDP.

[This work was supported in part by the Interdisciplinary Quantitative Biology (IQ Biology) program at the BioFrontiers Institute, University of Colorado, Boulder.]

## Repository contents

This repo contains the following scripts, described in the order with which they should be used.

### data_prep.R

#### Goal:

Prepare data for input into the disease model.

#### Inputs:

* **[2017_Gaz_zcta_national.txt](https://www.census.gov/geographies/reference-files/time-series/geo/gazetteer-files.2017.html#list-tab-264479560)** 
  * 2017 zip code latitude and longitudes from the US Census Bureau
* **ACSDP5Y2017.DP05-Data.csv** 
  * 2017 populations by zip code from the US Census Bureau
* **state_abbrev.csv** 
  * US states with their 2-letter abbreviations
* **migrants_merger.csv** 
  * demographic data for agricultural workforce
  * output from: [migrants_R_proj/migrants_stats.R at master · naiacasina/migrants_R_proj · GitHub](https://github.com/naiacasina/migrants_R_proj/blob/master/migrants_stats.R)

#### Outputs:

* **2017_pop_lat_long_data_states.csv** 
  * 3-digit zip: populations, lat, long, and state abbreviations
* **2017_pop_demo_data_agricultural.csv** 
  * 3-digit zip: populations, lat, long, state abbreviations, and agricultural worker demographic information



### parameterization.R

#### Goal:

Determine optimal values for the coefficients for demographic variables for input into the agricultural worker SIR model.

#### Inputs:

* **general_population_demographics.csv**
  * 2017 information on [household crowding](https://data.census.gov/table/ACSDP1Y2017.DP04?g=010XX00US$0400000) and [households with at least one child](https://data.census.gov/table/ACSST1Y2017.S1101?t=Families%20and%20Living%20Arrangements&g=010XX00US$0400000) from the US Census Bureau

#### Outputs:

* Optimal parameter values (printed to screen)



### R0_exploration.R

#### Goal:

Check to make sure the R0 values for the agricultural workforce would be within a reasonable range using the optimized parameter values from parameterization.R.

#### Inputs:

* optimized parameter values from parameterization.R

#### Outputs:

* summary statistics and plots for R0 values for the agricultural workforce based on optimized parameters (printed to screen)



### disease_model.R

#### Goal:

Develop and run the disease model and output model results.

#### Inputs:

* **2017_pop_demo_data_agricultural.csv**
  * 3-digit zip: populations, lat, long, state abbreviations, and agricultural worker demographic information
  * output from data_prep.R

#### Outputs:

* **infection_start_data.csv**
  * infection start dates for each 3-digit zip location
* **model_output_state_demographics.csv**
  * SIR model data with corresponding state information
* **model_output_grouped_by_state.csv**
  * SIR model data grouped by state for each timestep



### figures.R

#### Goal:

Plot disease model output and agricultural workforce data.

#### Inputs:

* **infection_start_data.csv**
  * infection start dates for each 3-digit zip location
  * output from disease_model.R
* **model_output_state_demographics.csv**
  * SIR model data with corresponding state information
  * output from disease_model.R
* **model_output_grouped_by_state.csv**
  * SIR model data grouped by state for each timestep
  * output from disease_model.R
* **migrants_merger.csv**
  * demographic data for agricultural workforce
  * output from: [migrants_R_proj/migrants_stats.R at master · naiacasina/migrants_R_proj · GitHub](https://github.com/naiacasina/migrants_R_proj/blob/master/migrants_stats.R)
* **general_population_demographics.csv**
  * 2017 information on [household crowding](https://data.census.gov/table/ACSDP1Y2017.DP04?g=010XX00US$0400000) and [households with at least one child](https://data.census.gov/table/ACSST1Y2017.S1101?t=Families%20and%20Living%20Arrangements&g=010XX00US$0400000) from the US Census Bureau

#### Outputs:

* pngs for the following figures saved to a "figures" directory:
  * map of spatial infection spread (infection_map.png)
  * histograms of demographic data (prop_crowded_hist.png; prop_children_hist.png)
  * maps of agricultural worker demographics (crowding_map.png; children_map.png; ag_workers_map.png)
  * peak infectious by state, colored by region (peak_infection_state_bar.png)
  * SIR plot of infectious by state, colored by region (infectious_by_state.png)

