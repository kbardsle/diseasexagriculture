ag_demographic_data <- read_csv("data/migrants_merger.csv")

ag_clean_demo_data <- demographic_data %>%
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


# read in demographic data
gen_demographic_data<- read_csv("data/general_population_demographics.csv")

# normalize demographic data
gen_demographic_data_normalized <- gen_demographic_data %>%
  mutate(proportion_crowded = percent_crowded/100,
         proportion_w_kids = percent_with_children/100)

ggplot(mapping=aes(x = proportion_crowded)) +
  geom_histogram(data = ag_clean_demo_data, fill="red", alpha=0.5) +
  geom_histogram(data = gen_demographic_data_normalized, fill="blue", alpha=0.5)

ggplot(mapping=aes(x = proportion_w_kids)) +
  geom_histogram(data = ag_clean_demo_data, fill="red", alpha=0.5) +
  geom_histogram(data = gen_demographic_data_normalized, fill="blue", alpha=0.5)
