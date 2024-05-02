ag_demographic_data <- read_csv("data/migrants_merger.csv")

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


# read in demographic data
gen_demographic_data<- read_csv("data/general_population_demographics.csv")

# normalize demographic data
gen_demographic_data_normalized <- gen_demographic_data %>%
  mutate(proportion_crowded = percent_crowded/100,
         proportion_w_kids = percent_with_children/100)

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
