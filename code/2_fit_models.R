library(tidyverse)
library(janitor)
library(brms)


# load data
thule_birds_clean = readRDS(file = "data/thule_birds_clean.rds") 

thule_birds_clean %>% 
  ggplot(aes(x = julian_z, y = hg_mean_species, color = species_common)) + 
  geom_point()

# linear regression - simple start

brm_lm = brm(hg ~ julian,
             data = thule_birds_clean,
             family = gaussian())


plot(conditional_effects(brm_lm), points = T)

# linear regression - z-score y axis
brm_lm_z = update(brm_lm, formula = . ~ julian_z,
                  newdata = thule_birds_clean)

plot(conditional_effects(brm_lm_z), points = T)


# linear regression - z-score y axis + species
brm_lm_z_species = update(brm_lm, formula = . ~ julian_z*species_common,
                  newdata = thule_birds_clean)

plot(conditional_effects(brm_lm_z_species), points = T)


# linear regression - gamma (positive real numbers. NO negative numbers. Higher mean = higher variance)

brm_lm_z_species = update(brm_lm, formula = . ~ julian_z*species_common,
                          newdata = thule_birds_clean,
                          family = Gamma(link = "log"))


plot(conditional_effects(brm_lm_z_species), points = T)


# generalized additive model with varying intercepts by species
# informative priors for intercept of normal(1, 1)

brm_smooth = brm(hg_mean_species ~ s(julian_z, by = species_common) + (1|species_common),
                 data = thule_birds_clean,
                 family = Gamma(link = "log"),
                 prior = c(prior(normal(0, 1), class = "Intercept"),
                           prior(exponential(2), class = "sd")),
                 iter = 2000,
                 chains = 4)


saveRDS(brm_smooth, file = "models/brm_smooth.rds")

plot_brm = plot(conditional_effects(brm_smooth_n_mi), points = T)

plot_brm$`hgmeanspecies.hgmeanspecies_julian_z:species_common` +
  facet_wrap(~species_common)

brm_smooth_n = update(brm_smooth, formula = . ~ s(julian_z, by = species_common) + 
                        (1|species_common) + nitrogen_s,
                      newdata = thule_birds_clean,
                      iter = 1000, chains = 1)
saveRDS(brm_smooth_n, file = "models/brm_smooth_n.rds")

plot(conditional_effects(brm_smooth_n), points = T)

brm_smooth_no_n = update(brm_smooth_n, formula = . ~ s(julian_z, by = species_common) + 
                           (1|species_common),
                         newdata = thule_birds_clean)


saveRDS(brm_smooth_no_n, file = "models/brm_smooth_no_n.rds")
plot(conditional_effects(brm_smooth_no_n), points = T)

waic(brm_smooth_n)
waic(brm_smooth_no_n)
waic(brm_smooth)

get_prior(bf(hg_mean_species ~ s(julian_z, by = species_common) + 
               (1|species_common) + mi(nitrogen_s)) + 
            bf(nitrogen_s|mi() ~ species_common),
          data = thule_birds_clean,
          prior = c(prior(normal(0, 1), class = "Intercept"),
                    prior(exponential(2), class = sd)))

brm_smooth_n_mi = brm(bf(hg_mean_species ~ s(julian_z, by = species_common) + 
                        (1|species_common) + mi(nitrogen_s)) + 
                          bf(nitrogen_s|mi() ~ species_common),
                      data = thule_birds_clean,
                      prior = c(prior(normal(0, 1), class = "Intercept")),
                      iter = 2000, chains = 4)

saveRDS(brm_smooth_n_mi, file = "models/brm_smooth_n_mi.rds")

brm_smooth_nc_mi = brm(bf(hg_mean_species ~ s(julian_z, by = species_common) + 
                           (1|species_common) + mi(nitrogen_s) + mi(carbon_s)) + 
                        bf(nitrogen_s|mi() ~ species_common) + 
                         bf(carbon_s|mi() ~ species_common),
                      data = thule_birds_clean,
                      prior = c(prior(normal(0, 1), class = "Intercept")),
                      iter = 1000, chains = 1)

saveRDS(brm_smooth_nc_mi, file = "models/brm_smooth_nc_mi.rds")

brm_smooth_n_wt = brm(bf(hg_mean_species ~ s(julian_z, by = species_common) + 
                            (1|species_common) + nitrogen_s + actual_weight_s),
                       data = thule_birds_clean,
                       prior = c(prior(normal(0, 1), class = "Intercept")),
                       iter = 1000, chains = 1)

saveRDS(brm_smooth_n_wt, file = "models/brm_smooth_n_wt.rds")


# compare models ----------------------------------------------------------

brm_smooth_nc_mi = readRDS(file = "models/brm_smooth_nc_mi.rds")
brm_smooth_n_mi = readRDS(file = "models/brm_smooth_n_mi.rds")
brm_smooth = readRDS(file = "models/brm_smooth.rds")
brm_smooth_no_n = readRDS(file = "models/brm_smooth_no_n.rds")
brm_smooth_n = readRDS(file = "models/brm_smooth_n.rds")
brm_smooth_n_wt = readRDS(file = "models/brm_smooth_n_wt.rds")


waic(brm_smooth_n_mi, newdata = thule_birds_clean %>% filter(!is.na(nitrogen))%>% filter(!is.na(carbon)))
waic(brm_smooth_nc_mi, newdata = thule_birds_clean %>% filter(!is.na(nitrogen))%>% filter(!is.na(carbon)))
waic(brm_smooth, newdata = thule_birds_clean %>% filter(!is.na(nitrogen)) %>% filter(!is.na(carbon)))
waic(brm_smooth_n, newdata = thule_birds_clean %>% filter(!is.na(nitrogen)) %>% filter(!is.na(carbon)))
waic(brm_smooth_no_n, newdata = thule_birds_clean %>% filter(!is.na(nitrogen)) %>% filter(!is.na(carbon)))
waic(brm_smooth_n_wt, newdata = thule_birds_clean %>% filter(!is.na(nitrogen)) %>% filter(!is.na(carbon)))


pp_check(brm_smooth_n, stat = "median", type = "stat")
