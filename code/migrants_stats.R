# from: https://github.com/naiacasina/migrants_R_proj/blob/master/migrants_stats.R
# written by Naia Ormaza Zulueta

# In this script we:
#   1. Compute distributions of variables across states for the NAWS survey
#      using weights
#   2. Compute the merger of data SURVEY + CENSUS

rm(list=ls()) 
packages <- c("dplyr", "survey", "ggplot2", "tidyr", "usethis", "arrow")
lapply(packages, require, character=TRUE)
# for git
#setwd('~/Library/CloudStorage/OneDrive-UCB-O365/Projects/Ag Workforce/Migrants/migrants_R_proj')

# change path to 'Migrants' folder
#setwd('/Users/naiacasina/Library/CloudStorage/OneDrive-UCB-O365/Projects/Ag Workforce/Migrants')
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

naws_ae <- read.csv('Data/NAWS_A2E197.csv')
naws_fz <- read.csv('Data/NAWS_F2Y197.csv')

# merger
naws_full <- merge(naws_ae, naws_fz, by = "FWID")

# drop years and non-migrants
naws <- naws_full %>%
  filter(FY.x %in% c(2022, 2017, 2012))
  # if we wnat to filter migrants, uncomment next row
  # %>% filter(MIGRANT %in% c(1))

# drop cols
threshold <- nrow(naws) / 2

# cleaned df
naws_df <- naws %>%
  select_if(~sum(is.na(.)) < threshold)

naws_df <- mutate(naws_df, PWTYCRD = ifelse(is.na(PWTYCRD.x), PWTYCRD.y, PWTYCRD.x))
naws_design <- svydesign(ids = ~1, data = naws_df, weights = ~PWTYCRD)

average_age <- svyby(~AGE, ~REGION6, naws_design, svymean)
gender_distribution <- svyby(~GENDER, ~REGION6, naws_design, svymean)

ggplot(gender_distribution, aes(x = REGION6, y = GENDER, fill = REGION6)) +
  geom_bar(stat = "identity") +
  labs(x = "Region", y = "Proportion of Females", title = "Gender Distribution by Region") +
  theme_minimal()

b17_region_distribution <- svytable(~B17CODE + REGION6, naws_design)

# test plot
b17_region_df <- as.data.frame(b17_region_distribution)
b17_region_df$Proportion <- b17_region_df$Freq / sum(b17_region_df$Freq)
b17_region_df <- b17_region_df %>%
  group_by(REGION6) %>%
  mutate(Region_Proportion = Freq / sum(Freq)) %>%
  ungroup()

ggplot(b17_region_df, aes(x = B17CODE, y = Region_Proportion, fill = B17CODE)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~REGION6) +  # Facet by REGION6 to separate plots by region
  labs(x = "Country of Origin", y = "Proportion within Region", 
       title = "Distribution of Country of Origin by Survey Region") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# -----------
# --------- Building a dataframe -------
# -----------

compute_global_bins <- function(data, variable) {
  # ranges across all data
  min_val <- min(data[[variable]], na.rm = TRUE)
  max_val <- max(data[[variable]], na.rm = TRUE)
  
  # define breaks for binning - creating 5 bins but can be modified
  breaks <- seq(from = min_val, to = max_val, length.out = 6) 
  return(breaks)
}

compute_proportions <- function(variable, data, weights_col, regions, year_col, breaks, continuous_vars, migrant_col) {
  results_df <- data.frame()
  
  for (region in regions) {
    for (year in unique(data[[year_col]])) {
      for (migrant_status in unique(data[[migrant_col]])) {
        subset_data <- data[data$REGION6 == region & data[[year_col]] == year & data[[migrant_col]] == migrant_status & !is.na(data[[variable]]), ]
        
        if (variable %in% continuous_vars) {
          subset_data$Category <- cut(subset_data[[variable]], breaks = breaks[[variable]], include.lowest = TRUE, right = FALSE)
        } else {
          subset_data$Category <- as.factor(subset_data[[variable]])
        }
        
        category_weights <- tapply(subset_data[[weights_col]], subset_data$Category, sum)
        total_weight <- sum(category_weights, na.rm = TRUE)
        proportions <- category_weights / total_weight
        
        # category data
        for (cat in names(proportions)) {
          results_df <- rbind(results_df, data.frame(
            REGION6 = region,
            FY = year,
            MIGRANT = migrant_status,
            Variable = paste(variable, gsub("\\((.*),.*", "\\1", cat), sep = "."),
            Proportion = proportions[cat]
          ))
        }
      }
    }
  }
  
  return(results_df)
}

categorical_vars <- c("B14CODE","B17CODE", "A21Cx", "A22Cx", "A24a", "A24b",
                      "HHFAMGRD", "HHGRDKID", "HHKID", "HHOTHFAM",
                      "HHPARENT", "HHSIB", "HHYTH018", "K018USFW",
                      "K018USNF","CROP", "D22", "D23", "D26", "D30", "D33A",
                      "D34Ax", "D35x", "D37A", "E01x", "G01")

binary_vars <- c("FTC", "GENDER", "INDIGENOUS", "INTLSHTL", "MARRIED",
                 "FAMPOV", "SPOUSE", 'SPOUSEFW', 'SPOUSENF', 'YOUTH',
                 'BLWAGE', "CROWDED1","CROWDED2", "D11", "D12WG4")

all_vars <- c(binary_vars, categorical_vars)

continuous_vars <- c("C09WEEKS", "FWRDAYS", "FWWEEKS")

# Compute bins for each continuous variable
global_breaks <- list()
for (var in continuous_vars) {
  global_breaks[[var]] <- compute_global_bins(naws_df, var)
}

year_col = "FY.x"
migrant_col = "MIGRANT"
all_results <- data.frame()

for (var in continuous_vars) {
  var_results <- compute_proportions(var, naws_df, "PWTYCRD.x", sort(unique(naws_df$REGION6)), "FY.x", global_breaks, continuous_vars, migrant_col)
  all_results <- rbind(all_results, var_results)
}

for (var in all_vars) {
  var_results <- compute_proportions(var, naws_df, "PWTYCRD.x", unique(naws_df$REGION6), "FY.x", global_breaks, continuous_vars, migrant_col)
  all_results <- rbind(all_results, var_results)
}

all_results <- all_results[!is.na(all_results$MIGRANT), ]

# convert to wide format
all_results_wide <- pivot_wider(all_results, names_from = Variable, values_from = Proportion, id_cols = c("REGION6", "FY", "MIGRANT"))

# include states to regions
region_mapping <- data.frame(
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
)

results_with_regions <- merge(all_results_wide, region_mapping, by = "REGION6", all.x = TRUE)

write_parquet(results_with_regions, 'Results/results_with_regions.parquet')

# --------------- merger ------------------
census_df <- read.csv('Data/61D69F31-5A33-3017-A901-3F34975714EB.csv')
census_df <- census_df %>%
  filter(Year %in% c(2017, 2012))
results_with_regions$States <- strsplit(as.character(results_with_regions$States), ",\\s*")

# create a new dataframe where each state has its own row with the corresponding REGION6 data
results_expanded <- results_with_regions %>%
  tidyr::unnest(States) %>%
  rename(State = States) %>%
  distinct()

# trim whitespace from the State names
results_expanded$State <- toupper(results_expanded$State)

census_df$Year <- as.integer(census_df$Year)
results_expanded$FY <- as.integer(results_expanded$FY)

census_df <- census_df %>%
  mutate(Category = case_when(
    grepl("HIRED", Data.Item) | grepl("UNPAID", Data.Item) ~ "Non-Migrant",
    grepl("MIGRANT", Data.Item) ~ "Migrant",
    TRUE ~ as.character(NA)
  )) %>%
  filter(!is.na(Category))

results_expanded <- results_expanded %>%
  mutate(Category = factor(MIGRANT, labels = c("Non-Migrant", "Migrant")),
         FY = as.character(FY))

# Rename columns for consistency
census_df <- rename(census_df, Year = Year, State = State)

merged_data <- merge(results_expanded, census_df, by.x = c("State", "FY", "Category"), by.y = c("State", "Year", "Category"), all = TRUE)

cleaned_data <- merged_data %>%
  select(where(~ !all(is.na(.))))

cleaned_data <- cleaned_data %>%
  filter(State != "ALASKA" & State != "HAWAII")

wide_data <- cleaned_data %>%
  pivot_wider(
    names_from = Domain.Category, 
    values_from = c("Value", "CV...."),
    names_sep = "_"
  )

wide_data <- wide_data %>%
  select(-c(MIGRANT, Domain, `Data.Item`, Commodity, watershed_code, `State.ANSI`, `Geo.Level`, Program, Period)) %>%
  rename(Value = 'Value_NOT SPECIFIED') %>%
  relocate(Value, .after = 3)  # Move "Value" to be the 3rd column


write_parquet(wide_data, 'Results/migrants_merger.parquet')
write.csv(wide_data, 'Data/migrants_merger.csv')

