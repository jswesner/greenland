library(tidyverse)
library(janitor)
library(brms)

# load data
thule_birds_clean = readRDS(file = "data/thule_birds_clean.rds") %>% 
  mutate(species_abbreviation = case_when(species_common == "Dovekie" ~ "e) DOVE",
                                          species_common == "Black-legged Kittiwake" ~ "d) BLKI",
                                          species_common == "Atlantic Puffin" ~ "c) ATPU",
                                          species_common == "Black Guillemot" ~ "b) BLGU",
                                          species_common == "Thick-billed Murre" ~ "a) TBMU"))

# load models
brm_smooth_nc_gamma = readRDS(file = "models/brm_smooth_nc_gamma.rds")
brm_smooth_nc_gamma_wide = readRDS(file = "models/brm_smooth_nc_gamma_wide.rds")

# brm_smooth_nc_gamma_wide = update(brm_smooth_nc_gamma,
#                                   prior = c(prior(normal(0, 1), class = "Intercept"),
#                                             prior(exponential(5), class = "sd"),
#                                             prior(normal(0, 1), class = "b"),
#                                             prior(normal(0, 2), class = "b", coef = "nitrogen_s")))
# 
# saveRDS(brm_smooth_nc_gamma_wide, file = "models/brm_smooth_nc_gamma_wide.rds")


# compare priors for time series ------------------------------------------

pred_hg_grid = thule_birds_clean %>% 
  group_by(species_common, sd_julian, mean_julian) %>% 
  mutate(min = min(julian_z),
         max = max(julian_z)) %>% 
  select(-julian_z) %>% 
  distinct(species_common, sd_julian, mean_julian, min, max, mean_hg) %>% 
  group_by(species_common, sd_julian, mean_julian, min, max, mean_hg) %>% 
  reframe(julian_z = seq(min, max, length.out = 30)) %>% 
  expand_grid(nitrogen_s = c(quantile(thule_birds_clean$nitrogen_s, na.rm = T))) %>% 
  expand_grid(carbon_s = c(quantile(thule_birds_clean$carbon_s, na.rm = T))) %>% 
  mutate(julian = (julian_z*sd_julian) + mean_julian,
         date = as.Date(julian)) %>% 
  mutate(species_abbreviation = case_when(species_common == "Dovekie" ~ "e) DOVE",
                                          species_common == "Black-legged Kittiwake" ~ "d) BLKI",
                                          species_common == "Atlantic Puffin" ~ "c) ATPU",
                                          species_common == "Black Guillemot" ~ "b) BLGU",
                                          species_common == "Thick-billed Murre" ~ "a) TBMU"))

post_preds_hg_original = pred_hg_grid %>% 
  add_epred_draws(brm_smooth_nc_gamma, re_formula = NULL) %>% 
  mutate(hg = .epred*mean_hg) %>% 
  group_by(species_common, species_abbreviation, .draw, julian, date, julian_z, mean_hg) %>% 
  reframe(hg = mean(hg),
          .epred = mean(.epred)) %>% 
  mutate(model_run = "original")

post_preds_hg_wide = pred_hg_grid %>% 
  add_epred_draws(brm_smooth_nc_gamma_wide, re_formula = NULL) %>% 
  mutate(hg = .epred*mean_hg) %>% 
  group_by(species_common, species_abbreviation, .draw, julian, date, julian_z, mean_hg) %>% 
  reframe(hg = mean(hg),
          .epred = mean(.epred)) %>% 
  mutate(model_run = "wide")

post_preds_hg_original_wide = bind_rows(post_preds_hg_original,
                                        post_preds_hg_wide)

post_timeseries_original_wide = post_preds_hg_original_wide %>% 
  ggplot(aes(x = date, y = hg)) + 
  stat_lineribbon(.width = c(0.5, 0.75, 0.95), alpha = 0.2, linewidth = 0.1,
                  aes(fill = species_abbreviation)) + 
  ggh4x::facet_grid2(model_run~species_abbreviation) +
  guides(fill = "none",
         color = "none") +
  theme_default() +
  labs(y = "Hg (ng/g wet blood)",
       x = "Date") +
  theme(strip.text = element_text(size = 8, hjust = 0),
        text = element_text(size = 8)) +
  scale_fill_colorblind()

ggsave(post_timeseries_original_wide, file = "plots/post_timeseries_original_wide.jpg",
       width = 6.5, height = 4)

# compare prior for nitrogen regression -----------------------------------
post_n_species_original = thule_birds_clean %>%
  select(species_common, nitrogen_s, sd_n, mean_n, mean_hg) %>% 
  group_by(species_common, sd_n, mean_n, mean_hg) %>% 
  mutate(min = min(nitrogen_s, na.rm = T),
         max = max(nitrogen_s, na.rm = T)) %>% 
  distinct(species_common, min, max, sd_n, mean_n, mean_hg) %>% 
  reframe(nitrogen_s = seq(min, max, length.out = 30)) %>% 
  mutate(species_abbreviation = case_when(species_common == "Dovekie" ~ "e) DOVE",
                                          species_common == "Black-legged Kittiwake" ~ "d) BLKI",
                                          species_common == "Atlantic Puffin" ~ "c) ATPU",
                                          species_common == "Black Guillemot" ~ "b) BLGU",
                                          species_common == "Thick-billed Murre" ~ "a) TBMU")) %>%
  expand_grid(julian_z = quantile(thule_birds_clean$julian_z, na.rm = T)) %>% 
  expand_grid(carbon_s = quantile(thule_birds_clean$carbon_s, na.rm = T))  %>% 
  mutate(nitrogen = (nitrogen_s*sd_n) + mean_n) %>% 
  add_epred_draws(brm_smooth_nc_gamma, 
                  re_formula = ~ (1 + nitrogen_s + carbon_s|species_common)) %>% 
  mutate(hg = .epred*mean_hg) %>% 
  group_by(.draw, nitrogen, species_common, species_abbreviation) %>% 
  reframe(hg = mean(hg),
          .epred = mean(.epred))

post_n_species_wide = thule_birds_clean %>%
  select(species_common, nitrogen_s, sd_n, mean_n, mean_hg) %>% 
  group_by(species_common, sd_n, mean_n, mean_hg) %>% 
  mutate(min = min(nitrogen_s, na.rm = T),
         max = max(nitrogen_s, na.rm = T)) %>% 
  distinct(species_common, min, max, sd_n, mean_n, mean_hg) %>% 
  reframe(nitrogen_s = seq(min, max, length.out = 30)) %>% 
  mutate(species_abbreviation = case_when(species_common == "Dovekie" ~ "e) DOVE",
                                          species_common == "Black-legged Kittiwake" ~ "d) BLKI",
                                          species_common == "Atlantic Puffin" ~ "c) ATPU",
                                          species_common == "Black Guillemot" ~ "b) BLGU",
                                          species_common == "Thick-billed Murre" ~ "a) TBMU")) %>%
  expand_grid(julian_z = quantile(thule_birds_clean$julian_z, na.rm = T)) %>% 
  expand_grid(carbon_s = quantile(thule_birds_clean$carbon_s, na.rm = T))  %>% 
  mutate(nitrogen = (nitrogen_s*sd_n) + mean_n) %>% 
  add_epred_draws(brm_smooth_nc_gamma_wide, 
                  re_formula = ~ (1 + nitrogen_s + carbon_s|species_common)) %>% 
  mutate(hg = .epred*mean_hg) %>% 
  group_by(.draw, nitrogen, species_common, species_abbreviation) %>% 
  reframe(hg = mean(hg),
          .epred = mean(.epred))


post_original_wide = bind_rows(post_n_species_original %>% mutate(model_run = "original priors"),
          post_n_species_wide %>% mutate(model_run = "wide priors"))

plot_n_original_wide = post_original_wide %>% 
  ggplot(aes(x = nitrogen, y = hg)) +
  stat_lineribbon(.width = c(0.5, 0.75, 0.95), aes(fill = species_abbreviation),
                  alpha = 0.2, linewidth = 0.2) +
  ggh4x::facet_grid2(model_run~species_abbreviation) +
  labs(y = "Hg (ng/g wet blood)",
       x = expression(delta^15 * N)) +
  scale_fill_colorblind() +
  guides(fill = "none") +
  theme(strip.text = element_text(size = 8, hjust = 0),
        text = element_text(size = 8)) +
  NULL

ggsave(plot_n_original_wide, file = "plots/plot_n_original_wide.jpg",
       width = 6.5, height = 4)
