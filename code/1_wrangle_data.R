library(tidyverse)
library(janitor)

thule_birds = read_csv("data/thule_birds_r.csv") %>% 
  clean_names() %>% 
  mutate(date = dmy(date),
         julian = julian(date)) %>% 
  mutate(mean_julian = mean(julian),
         sd_julian = sd(julian),
         julian_z = (julian - mean_julian)/sd_julian,
         hg_mean = hg/mean(hg)) %>% 
  group_by(species_common) %>% 
  mutate(hg_mean_species = hg/mean(hg),
         mean_n = mean(nitrogen, na.rm = T),
         sd_n = sd(nitrogen, na.rm = T),
         nitrogen_s = (nitrogen - mean_n)/sd_n,
         mean_c = mean(carbon, na.rm = T),
         sd_c = sd(carbon, na.rm = T),
         carbon_s = (carbon - mean_c)/sd_c,
         mean_hg = mean(hg)) %>% 
  mutate(mean_julian = mean(julian),
         sd_julian = sd(julian),
         julian_z = (julian - mean_julian)/sd_julian,
         hg_mean = hg/mean(hg),
         actual_weight = parse_number(actual_weight),
         actual_weight_s = scale(actual_weight)) 

saveRDS(thule_birds, file = "data/thule_birds_clean.rds")


thule_birds %>% 
  ggplot(aes(x = carbon, y = nitrogen, color = species_common)) +
  geom_point() +
  # geom_smooth() +
  # facet_wrap(~species_common) +
  NULL

thule_birds %>% 
  ggplot(aes(x = nitrogen_s, y = actual_weight)) + 
  geom_point() +
  facet_wrap(~species_common, scales = "free_x")
