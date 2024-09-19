library(tidyverse)
library(tidybayes)
library(brms)

# load models
brm_smooth_nc_gamma = readRDS("models/brm_smooth_nc_gamma.rds")
brm_smooth_isotope_model = readRDS(file = "models/brm_smooth_isotope_model.rds")

# make data grid to predict over
species_means = thule_birds_clean %>% distinct(species_common, mean_hg)
mean_julian = unique(thule_birds_clean$mean_julian)
sd_julian = unique(thule_birds_clean$sd_julian)

get_group_seqs = function(data, group, seq_var, length_out){
  temp = data %>% 
    group_by(group) %>% 
    group_split()
  
  group_seq = seq(min(temp$seq_var),
                  max(temp$seq_var),
                  length.out = length_out)
}

group_seqs = NULL
temp = thule_birds_clean %>% group_by(species_common) %>% group_split()

for(i in 1:length(temp)){
  group_seqs[[i]] = tibble(julian_z = seq(min(temp[[i]]$julian_z, na.rm = T),
                                          max(temp[[i]]$julian_z, na.rm = T),
                                          length.out = 30),
                           nitrogen_s = seq(min(temp[[i]]$nitrogen_s, na.rm = T),
                                            max(temp[[i]]$nitrogen_s, na.rm = T),
                                            length.out = 30)) %>% 
    mutate(species_common = unique(temp[[i]]$species_common))
  
}

# plot time series --------------------------------------------------------

post_preds_hg = bind_rows(group_seqs) %>%
  select(-nitrogen_s) %>% 
  expand_grid(nitrogen_s = c(quantile(thule_birds_clean$nitrogen_s, na.rm = T))) %>% 
  expand_grid(carbon_s = c(quantile(thule_birds_clean$carbon_s, na.rm = T))) %>% 
  mutate(julian = (julian_z*sd_julian) + mean_julian,
         date = as.Date(julian)) %>% 
  mutate(species_abbreviation = case_when(species_common == "Dovekie" ~ "e) DOVE",
                                          species_common == "Black-legged Kittiwake" ~ "d) BLKI",
                                          species_common == "Atlantic Puffin" ~ "c) ATPU",
                                          species_common == "Black Guillemot" ~ "b) BLGU",
                                          species_common == "Thick-billed Murre" ~ "a) TBMU")) %>% 
  left_join(species_means) %>% 
  add_epred_draws(brm_smooth_nc_gamma, re_formula = NULL) %>% 
  mutate(hg = .epred*mean_hg) %>% 
  group_by(species_common, species_abbreviation, .draw, julian, date, julian_z, mean_hg) %>% 
  reframe(hg = mean(hg),
          .epred = mean(.epred))

saveRDS(post_preds_hg, file = "posteriors/post_preds_hg.rds")
saveRDS(group_seqs, file = "data/group_seqs.rds")



# isotope time series data ------------------------------------------------
isotope_data = readRDS("data/isotope_data.rds") %>% ungroup

post_isotopes = isotope_data %>% 
  group_by(species_isotope, isotope_mean) %>% 
  reframe(julian_z = seq(min(julian_z), 
                         max(julian_z),
                         length.out = 30)) %>% 
  left_join(isotope_data %>% ungroup %>% distinct(species_isotope, mean_julian, sd_julian)) %>% 
  mutate(julian = (julian_z*sd_julian) + mean_julian,
         date = as.Date(julian)) %>% 
  separate(species_isotope, into = c("species_common", "isotope", NA),
           sep = "_", remove = F) %>% 
  add_epred_draws(brm_smooth_isotope_model) %>%  
  mutate(.epred = .epred*isotope_mean) 
      
saveRDS(post_isotopes, file = "posteriors/post_isotopes.rds")
