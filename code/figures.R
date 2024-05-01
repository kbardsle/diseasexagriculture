
# read in data (alternatively, get from disease model output directly)

fig_data <- read.csv(output_state_info_df, "data/model_output_state_demographics.csv")
infection_data <- read.csv("data/infection_start_data.csv")



# PLOTTING -------------------------------------

pops_50 <- sample(unique(data_reformat$Population), 50)
pops_4 <- sample(unique(data_reformat$Population), 4)

# figure with facets for populations
output_state_info_df %>% filter(Population %in% pops_4) %>% 
  ggplot(aes(x = step, y = n, col = Disease_Status, linetype=Community)) +
  geom_line() + facet_wrap(~Population) + theme_bw()

# figure with all pops plotted in the same figure, color by disease status
output_state_info_df %>% 
  ggplot(aes(x = step, y = n, col = Disease_Status, linetype=Community, group_by=Population)) +
  geom_line() + theme_bw()

 # map plot
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


# map for household crowding
map_df <- data_2017 %>% 
  group_by(State) %>% 
  summarize(crowded=mean(proportion_crowded),
            children=mean(proportion_w_kids),
            state_ag_pop=mean(state_ag_population),
            ag_prop=(median(state_ag_population)/sum(POP3)))
map_df$State <- tolower(map_df$State)
colnames(map_df)[1] <- "region"

full_map_df <- full_join(us_map, map_df, by="region")

map_crowding <- ggplot() + 
  geom_polygon(data = full_map_df,
               aes(x = long, y = lat, group = group, fill = crowded)) +
  geom_polygon(data = map_data("state"), 
               aes(x = long, y = lat, group = group), 
               color = "black", fill = NA) +
  scale_fill_viridis() + 
  theme_minimal() + 
  labs(fill="Proportion", title="Proportion of households with more than one resident per room")

map_children <- ggplot() + 
  geom_polygon(data = full_map_df,
               aes(x = long, y = lat, group = group, fill = children)) +
  geom_polygon(data = map_data("state"), 
               aes(x = long, y = lat, group = group), 
               color = "black", fill = NA) +
  scale_fill_viridis() + 
  theme_minimal() + 
  labs(fill="Proportion", title="Proportion of households with children")

map_ag_workers <- ggplot() + 
  geom_polygon(data = full_map_df,
               aes(x = long, y = lat, group = group, fill = ag_prop)) +
  geom_polygon(data = map_data("state"), 
               aes(x = long, y = lat, group = group), 
               color = "black", fill = NA) +
  scale_fill_viridis() + 
  theme_minimal() +
  labs(fill="Proportion", title="Proportion of state population that are agricultural workers")


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