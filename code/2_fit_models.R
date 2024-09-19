library(tidyverse)
library(janitor)
library(brms)


# load data
thule_birds_clean = readRDS(file = "data/thule_birds_clean.rds") 

# prior predictive
brm_smooth_prior = brm(hg_mean_species ~ s(julian_z, by = species_common) +
                         nitrogen_s + (1 + nitrogen_s|species_common),
                 data = thule_birds_clean,
                 family = Gamma(link = "log"),
                 prior = c(prior(normal(0, 0.1), class = "Intercept"),
                           prior(exponential(5), class = "sd"),
                           prior(normal(0, 0.01), class = "b"),
                           prior(normal(0, 1), class = "b", coef = "nitrogen_s")),
                 iter = 500,
                 chains = 1,
                 sample_prior = "only")


saveRDS(brm_smooth_prior, file = "models/brm_smooth_prior.rds")

# plot and compare to risk cutoffs of 1000 ng/g/500. Dividing by 500 rescales the references to the mean-scaled Hg in the model
prior_conds = plot(conditional_effects(brm_smooth_prior, effects = "julian_z:species_common", prob = 0.5))

prior_conds$`julian_z:species_common` + scale_y_log10() + geom_hline(yintercept = c(1000/500, 2000/500))

brm_smooth_n_gamma = update(readRDS("models/brm_smooth_prior.rds"), formula = . ~ 1 + s(julian_z, by = species_common) + 
                        (1 + nitrogen_s|species_common),
                      newdata = thule_birds_clean,
                      iter = 2000, chains = 4, sample_prior = "no",
                      cores = 4)

saveRDS(brm_smooth_n_gamma, file = "models/brm_smooth_n_gamma.rds")


brm_smooth_c_gamma = update(readRDS("models/brm_smooth_prior.rds"), formula = . ~ 1 + s(julian_z, by = species_common) + 
                              (1 + carbon_s|species_common),
                            newdata = thule_birds_clean,
                            iter = 2000, chains = 4, sample_prior = "no",
                            cores = 4)

saveRDS(brm_smooth_c_gamma, file = "models/brm_smooth_c_gamma.rds")


brm_smooth_nc_gamma = update(readRDS("models/brm_smooth_n_gamma.rds"), 
                             formula = . ~ 1 + s(julian_z, by = species_common) + 
                              (1 + nitrogen_s + carbon_s|species_common),
                            newdata = thule_birds_clean,
                            iter = 2000, chains = 4, sample_prior = "no",
                            cores = 4)

saveRDS(brm_smooth_nc_gamma, file = "models/brm_smooth_nc_gamma.rds")


brm_smooth_no_nc_gamma = update(brm_smooth_n_gamma, formula = . ~ s(julian_z, by = species_common) + 
                           (1|species_common),
                           cores =4,
                         newdata = thule_birds_clean)


saveRDS(brm_smooth_no_nc_gamma, file = "models/brm_smooth_no_nc_gamma.rds")

# compare models ----------------------------------------------------------
brm_smooth_c_gamma = readRDS(file = "models/brm_smooth_c_gamma.rds")
brm_smooth_n_gamma = readRDS(file = "models/brm_smooth_n_gamma.rds")
brm_smooth_nc_gamma = readRDS(file = "models/brm_smooth_nc_gamma.rds")
brm_smooth_no_nc_gamma = readRDS(file = "models/brm_smooth_no_nc_gamma.rds")

waic_newdata = brm_smooth_nc_gamma$data

waic(brm_smooth_c_gamma, newdata = waic_newdata)
waic(brm_smooth_n_gamma, newdata = waic_newdata)
waic(brm_smooth_nc_gamma, newdata = waic_newdata)
waic(brm_smooth_no_nc_gamma, newdata = waic_newdata)

pp_check(brm_smooth_n_gamma, stat = "median", type = "stat")
pp_check(brm_smooth_c_gamma, stat = "median", type = "stat")
pp_check(brm_smooth_nc_gamma, stat = "median", type = "stat")
pp_check(brm_smooth_no_nc_gamma, stat = "median", type = "stat")

pp_check(brm_smooth_n_gamma)
pp_check(brm_smooth_c_gamma)
pp_check(brm_smooth_nc_gamma)
pp_check(brm_smooth_no_nc_gamma)


# fit with a globally centered response instead of centered by species -----------------------------------------------------

newdat_1 = thule_birds_clean %>% 
  ungroup %>% 
  mutate(hg_1 = hg/mean(hg, na.rm = T))

brm_smooth_nc_gamma_1 = update(readRDS("models/brm_smooth_n_gamma.rds"), 
                             formula = hg_1 ~ s(julian_z, by = species_common) + (1 + nitrogen_s + carbon_s|species_common),
                             newdata = newdat_1,
                             family = Gamma(link = "log"),
                             prior = c(prior(normal(0, 0.5), class = "Intercept"),
                                       prior(exponential(5), class = "sd"),
                                       prior(normal(0, 0.01), class = "b"),
                                       prior(normal(0, 1), class = "b", coef = "nitrogen_s")),
                             iter = 2000, chains = 4, sample_prior = "no",
                             cores = 4)

saveRDS(brm_smooth_nc_gamma_1, file = "models/brm_smooth_nc_gamma_1.rds")

# fit nitrogen and carbon models ------------------------------------------
isotope_data  = thule_birds_clean %>% 
  pivot_longer(cols = c(nitrogen, carbon),
               names_to = "isotope", 
               values_to = "isotope_values") %>%
  group_by(species_common, isotope) %>% 
  mutate(isotope_mean = mean(isotope_values, na.rm = T)) %>% 
  mutate(isotope_values_mean_centered = isotope_values/isotope_mean) %>% 
  mutate(species_isotope = paste0(species_common, "_", isotope))

saveRDS(isotope_data, file = "data/isotope_data.rds")

brm_smooth_isotope_model = brm(isotope_values_mean_centered ~ s(julian_z, by = species_isotope),
                               data = isotope_data, 
                               family = gaussian(),
                               prior = c(prior(normal(1, 0.2), class = "Intercept"),
                                         prior(normal(0, 0.1), class = "b")),
                               chains = 4, iter = 2000)

saveRDS(brm_smooth_isotope_model, file = "models/brm_smooth_isotope_model.rds")
