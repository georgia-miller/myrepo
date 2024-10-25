
## Script to analyse and plot cfu count

# This script reads a CSV file with CFU counts, processes them through finding the 
# average CFU count across spots, corrects for dilution factor, normalises to the 
# immediately post-infection time-point and plots the data.

# This script reads CSV files with the following columns:
# * date            (of counting final plate)
# * number          (ID so each condition and timepoint has a separate ID for validation)
# * spot_number     (usually 1-3)
# * strain          (salmonella),
# * infection_time  (post-infection)
# * condition       (pre-infection treatment)
# * antibiotic      (post-infection antibiotic used)
# * dilution        (factor that cfu was counted at e.g. 4 = 10^-4)
# * cfu_count       (per spot)

# For a new CFU experiment, copy the template section into a new script and replace:
#   * GKMXX
#   * Date (in plot title)

# You could then copy it back as a new section as all data frames are saved with 
# an experiment identifier (GKMXX) so can be analysed and plotted in parallel if needed.

# If needed, also change:
#   * Strain numbers (e.g. ST001)
#   * Condition order (keep in mind this is saved under the same variable "order" for all experiments)

# Remember to change your working directory and file paths.

library(dplyr)
library(ggplot2)

setwd("/PATH/TO/YOUR/DIRECTORY/")


# Process GKM02 data ------------------------------------------------------
##Load and process data
GKM02_data <- read.csv("cfu_analysis_example_data.csv") %>%
  # find average cfu across the three spots, grouped by all other factors so appear in new tibble
  group_by(number, strain, infection_time, condition, antibiotic, dilution) %>% 
  summarise(cfu_count = mean(cfu_count)) %>% 
  # find number of bacteria per ml, took 20 ul (*50)
  mutate(cfu_per_ml = cfu_count * 50) %>% 
  # correct for dilution by * 10^dilution factor
  mutate(cfu_corrected = cfu_per_ml * 10^dilution )


##Normalise 18 h against 1 h
# extract 1 h baseline for ST001
GKM02_baseline_st001 <-
  GKM02_data %>% 
  filter(infection_time == "1", strain == "st001") %>% 
  pull(cfu_corrected)

# normalise 18 h against this 1 h baseline ST001
GKM02_st001 <-
  GKM02_data %>% 
  filter(strain == "st001") %>% 
  mutate(cfu_norm = cfu_corrected / GKM02_baseline_st001)


# extract 1 h baseline for ST016
GKM02_baseline_st016 <-
  GKM02_data %>% 
  filter(infection_time == "1", strain == "st016") %>% 
  pull(cfu_corrected)

# normalise 18 h against this 1 h baseline ST016
GKM02_st016 <-
  GKM02_data %>% 
  filter(strain == "st016") %>% 
  mutate(cfu_norm = cfu_corrected / GKM02_baseline_st016)


# merge data
GKM02_data_full <-   
  full_join(GKM02_st001, GKM02_st016) %>% 
  # make cfu_norm into percentage
  mutate(cfu_percentage = cfu_norm * 100)


# Plot GKM02 data ---------------------------------------------------------

# plot initially, check the baseline is at 100% as expected then can filter out
GKM02_data_full %>% 
  ggplot(mapping = aes(x = condition, y = cfu_percentage, colour = strain)) +
  geom_point(size = 4, alpha = 0.7) # transparency so can check both baselines are at 100%

# set desired order of conditions for plot
order <- c("none", "tet", "doxy", "tet_doxy")

# now plot with more specifications
GKM02_plot <-
  GKM02_data_full %>% 
  filter(infection_time != "1") %>% 
  # reorder columns as specified by order variable (ensure data points also swap, not just the labels)
  mutate(condition = factor(condition, levels = order) ) %>%  
  
  ggplot(mapping = aes(x = condition, y = cfu_percentage, colour = strain)) +
  geom_point(size = 4, position = position_dodge(width = 0.75)) +
  labs(title = "GKM02 20241018", x = "Pre-treatment", y = "CFU 18 h / 1 h (%)")
#theme_classic()            # if want classic theme, add + at end of line above and uncomment

# visualise plot
GKM02_plot

# save plot, do not overwrite previous plots
ggsave("GKM02_plot.png", plot = GKM02_plot, path = "ENTER/PATH")



# TEMPLATE ----------------------------------------------------------------

# Process GKMXX data ------------------------------------------------------
##Load and process data
GKMXX_data <- read.csv("cfu_analysis_GKMXX_data.csv") %>%
  # find average cfu across the three spots, grouped by all other factors so appear in new tibble
  group_by(number, strain, infection_time, condition, antibiotic, dilution) %>% 
  summarise(cfu_count = mean(cfu_count)) %>% 
  # find number of bacteria per ml, took 20 ul (*50)
  mutate(cfu_per_ml = cfu_count * 50) %>% 
  # correct for dilution by * 10^dilution factor
  mutate(cfu_corrected = cfu_per_ml * 10^dilution )


##Normalise 18 h against 1 h
# extract 1 h baseline for ST001
GKMXX_baseline_st001 <-
  GKMXX_data %>% 
  filter(infection_time == "1", strain == "st001") %>% 
  pull(cfu_corrected)

# normalise 18 h against this 1 h baseline ST001
GKMXX_st001 <-
  GKMXX_data %>% 
  filter(strain == "st001") %>% 
  mutate(cfu_norm = cfu_corrected / GKMXX_baseline_st001)


# extract 1 h baseline for ST016
GKMXX_baseline_st016 <-
  GKMXX_data %>% 
  filter(infection_time == "1", strain == "st016") %>% 
  pull(cfu_corrected)

# normalise 18 h against this 1 h baseline ST016
GKMXX_st016 <-
  GKMXX_data %>% 
  filter(strain == "st016") %>% 
  mutate(cfu_norm = cfu_corrected / GKMXX_baseline_st016)


# merge data
GKMXX_data_full <-   
  full_join(GKMXX_st001, GKMXX_st016) %>% 
  # make cfu_norm into percentage
  mutate(cfu_percentage = cfu_norm * 100)


# Plot GKMXX data ---------------------------------------------------------

# plot initially, check the baseline is at 100% as expected then can filter out
GKMXX_data_full %>% 
  ggplot(mapping = aes(x = condition, y = cfu_percentage, colour = strain)) +
  geom_point(size = 4, alpha = 0.7) # transparency so can check both baselines are at 100%

# set desired order of conditions for plot
order <- c("none", "tet", "doxy", "tet_doxy")

# now plot with more specifications
GKMXX_plot <-
  GKMXX_data_full %>% 
  filter(infection_time != "1") %>% 
  # reorder columns as specified by order variable (ensure data points also swap, not just the labels)
  mutate(condition = factor(condition, levels = order) ) %>%  
  
  ggplot(mapping = aes(x = condition, y = cfu_percentage, colour = strain)) +
  geom_point(size = 4, position = position_dodge(width = 0.75)) +
  labs(title = "GKMXX 2024XXXX", x = "Pre-treatment", y = "CFU 18 h / 1 h (%)")
#theme_classic()            # if want classic theme, add + at end of line above and uncomment

# visualise plot
GKMXX_plot

# save plot, do not overwrite previous plots
ggsave("GKMXX_plot.png", plot = GKMXX_plot, path = "ENTER/PATH")


