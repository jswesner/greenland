library(tidyverse)
library(tidybayes)
library(ggthemes)
library(brms)
theme_set(brms::theme_default())

# load data
thule_birds_clean = readRDS(file = "data/thule_birds_clean.rds")  %>% 
  mutate(species_abbreviation = case_when(species_common == "Dovekie" ~ "e) DOVE",
                                          species_common == "Black-legged Kittiwake" ~ "d) BLKI",
                                          species_common == "Atlantic Puffin" ~ "c) ATPU",
                                          species_common == "Black Guillemot" ~ "b) BLGU",
                                          species_common == "Thick-billed Murre" ~ "a) TBMU"))


# load model
brm_smooth_nc_gamma = readRDS("models/brm_smooth_nc_gamma.rds")

# load posteriors
post_preds_hg = readRDS(file = "posteriors/post_preds_hg.rds") 

# plot time series --------------------------------------------------------
post_timeseries = post_preds_hg %>% 
  ggplot(aes(x = date, y = hg)) + 
  stat_lineribbon(.width = c(0.5, 0.75, 0.95), alpha = 0.2, linewidth = 0.1,
                  aes(fill = species_abbreviation)) + 
  facet_wrap(~species_abbreviation, nrow = 1) +
  geom_point(data = thule_birds_clean, aes(y = hg), shape = ".") +
  guides(fill = "none",
         color = "none") +
  theme_default() +
  labs(y = "Hg (ng/g wet blood)",
       x = "Date") +
  theme(strip.text = element_text(size = 8, hjust = 0),
        text = element_text(size = 8)) +
  scale_fill_colorblind()
  
ggsave(post_timeseries, file = "plots/post_timeseries.jpg", width = 6.5, height = 2, dpi = 400)


# species Hg --------------------------------------------------------------
post_means = post_preds_hg %>% 
  group_by(species_common) %>% 
  mutate(median = median(hg)) %>% 
  group_by(species_common, .draw, median) %>% 
  reframe(hg = mean(hg))

post_mean_plot = post_means %>% 
  ggplot(aes(x = reorder(species_common, median),
             y = hg)) +
  geom_point(data = thule_birds_clean %>% left_join(post_means %>% ungroup %>%
                                                      distinct(species_common, median)),
             aes(y = hg),
             shape = "|",
             size = 0.4) +
  stat_halfeye(size = 0.1, linewidth = 0.5) +
  guides(fill = "none",
         color = "none") +
  theme_default() +
  labs(y = "Hg (ng/g wet blood)",
       x = "") +
  theme(strip.text = element_text(size = 11),
        text = element_text(size = 11)) +
  scale_fill_colorblind() + 
  scale_color_colorblind() +
  # scale_y_log10() +
  geom_hline(yintercept = c(200, 1000), linewidth = 0.1, linetype = "dotted") +
  coord_flip() 

ggsave(post_mean_plot, file = "plots/post_mean_plot.jpg", width = 5, height = 3, dpi = 400)


# plot nitrogen -------------------------------------------
# by species
post_n_species = thule_birds_clean %>%
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


# plot n
  
saveRDS(post_n_species, file = "posteriors/post_n_species.rds")

plot_n_overall = post_n_species %>% 
  ggplot(aes(x = nitrogen, y = hg)) +
  geom_point(data = thule_birds_clean, aes(y = hg),
             shape = ".") + 
  stat_lineribbon(.width = c(0.5, 0.75, 0.95), aes(fill = species_abbreviation),
                  alpha = 0.2, linewidth = 0.2) +
  facet_wrap(~species_abbreviation, nrow = 1) +
  labs(y = "Hg (ng/g wet blood)",
       x = expression(delta^15 * N)) +
  scale_fill_colorblind() +
  guides(fill = "none") +
  theme(strip.text = element_text(size = 8, hjust = 0),
        text = element_text(size = 8)) +
  NULL

ggsave(plot_n_overall, file = "plots/plot_n_overall.jpg", 
       dpi = 400, width = 6.5, height = 2)


# plot carbon -------------------------------------------------------------

# by species
post_c_species = thule_birds_clean %>% ungroup %>% 
  group_by(species_common) %>% 
  mutate(min = min(carbon_s, na.rm = T),
         max = max(carbon_s, na.rm = T)) %>% 
  distinct(species_common, min, max) %>% 
  group_by(species_common) %>% 
  reframe(carbon_s = seq(min, max, length.out = 30 )) %>% 
  expand_grid(julian_z = quantile(thule_birds_clean$julian_z, na.rm = T)) %>% 
  expand_grid(nitrogen_s = quantile(thule_birds_clean$nitrogen_s, na.rm = T)) %>%
  add_epred_draws(brm_smooth_nc_gamma, 
                  re_formula = ~ (1 + nitrogen_s + carbon_s|species_common)) %>% 
  left_join(species_means) %>% 
  mutate(species_abbreviation = case_when(species_common == "Dovekie" ~ "e) DOVE",
                                          species_common == "Black-legged Kittiwake" ~ "d) BLKI",
                                          species_common == "Atlantic Puffin" ~ "c) ATPU",
                                          species_common == "Black Guillemot" ~ "b) BLGU",
                                          species_common == "Thick-billed Murre" ~ "a) TBMU")) %>% 
  mutate(carbon = (carbon_s*sd_carbon) + mean_carbon) %>% 
  mutate(hg = .epred*mean_hg) %>% 
  group_by(carbon_s, species_common, species_abbreviation, .draw, carbon) %>% 
  reframe(hg = mean(hg))


# plot c
saveRDS(post_c_species, file = "posteriors/post_c_species.rds")

plot_c_overall = post_c_species %>% 
  ggplot(aes(x = carbon, y = hg)) +
  geom_point(data = thule_birds_clean, aes(y = hg),
             shape = ".") + 
  stat_lineribbon(.width = c(0.5, 0.75, 0.95), aes(fill = species_abbreviation),
                  alpha = 0.2, linewidth = 0.2) +
  facet_wrap(~species_abbreviation, nrow = 1) +
  labs(y = "Hg (ng/g wet blood)",
       x = expression(delta^13 * C)) +
  scale_fill_colorblind() +
  guides(fill = "none") +
  theme(strip.text = element_text(size = 8, hjust = 0),
        text = element_text(size = 8)) +
  NULL

ggsave(plot_c_overall, file = "plots/plot_c_overall.jpg", 
       dpi = 400, width = 6.5, height = 2)

library(patchwork)
plot_n_and_c = plot_n_overall / plot_c_overall + theme(strip.text = element_blank())

ggsave(plot_n_and_c, file = "plots/plot_n_and_c.jpg", 
       dpi = 400, width = 6.5, height = 4)

# isotope plots -----------------------------------------------------------

isotope_data = readRDS("data/isotope_data.rds") %>% ungroup %>% 
  mutate(species_abbreviation = case_when(species_common == "Dovekie" ~ "e) DOVE",
                                          species_common == "Black-legged Kittiwake" ~ "d) BLKI",
                                          species_common == "Atlantic Puffin" ~ "c) ATPU",
                                          species_common == "Black Guillemot" ~ "b) BLGU",
                                          species_common == "Thick-billed Murre" ~ "a) TBMU"))

post_isotopes = readRDS("posteriors/post_isotopes.rds") %>% 
  mutate(species_abbreviation = case_when(species_common == "Dovekie" ~ "e) DOVE",
                                          species_common == "Black-legged Kittiwake" ~ "d) BLKI",
                                          species_common == "Atlantic Puffin" ~ "c) ATPU",
                                          species_common == "Black Guillemot" ~ "b) BLGU",
                                          species_common == "Thick-billed Murre" ~ "a) TBMU"))

library(ggh4x)
nitrogen_isotope_plot = post_isotopes %>% 
  filter(.draw <= 500) %>% 
  filter(isotope == "nitrogen") %>% 
  ggplot(aes(x = date, y = .epred)) + 
  geom_point(data = isotope_data %>% 
               filter(isotope == "nitrogen") , aes(y = isotope_values),
             shape = ".") +
  stat_lineribbon(.width = c(0.5, 0.75, 0.95), alpha = 0.2, aes(fill = species_common),
                  linewidth = 0.2) +
  facet_grid2(~species_abbreviation, scales = "free_y") +
  guides(fill = "none",
         color = "none") +
  theme_default() +
  labs(y = expression(delta^15 * N),
       x = "Date") +
  theme(strip.text = element_text(size = 8, hjust = 0),
        text = element_text(size = 8)) +
  scale_fill_colorblind()

isotope_data_carbon = readRDS("data/isotope_data.rds") %>% ungroup %>% 
  filter(isotope == "carbon") %>% 
  mutate(species_abbreviation = case_when(species_common == "Dovekie" ~ "j)",
                                          species_common == "Black-legged Kittiwake" ~ "i)",
                                          species_common == "Atlantic Puffin" ~ "h)",
                                          species_common == "Black Guillemot" ~ "g)",
                                          species_common == "Thick-billed Murre" ~ "f)"))

post_isotopes_carbon = readRDS("posteriors/post_isotopes.rds") %>% 
  filter(isotope == "carbon") %>% 
  mutate(species_abbreviation = case_when(species_common == "Dovekie" ~ "j)",
                                          species_common == "Black-legged Kittiwake" ~ "i)",
                                          species_common == "Atlantic Puffin" ~ "h)",
                                          species_common == "Black Guillemot" ~ "g)",
                                          species_common == "Thick-billed Murre" ~ "f)"))

carbon_isotope_plot = post_isotopes_carbon %>% 
  filter(.draw <= 500) %>% 
  ggplot(aes(x = date, y = .epred)) + 
  geom_point(data = isotope_data_carbon , aes(y = isotope_values),
             shape = ".") +
  stat_lineribbon(.width = c(0.5, 0.75, 0.95), alpha = 0.2, aes(fill = species_common),
                  linewidth = 0.2) +
  facet_grid2(~species_abbreviation, scales = "free_y") +
  guides(fill = "none",
         color = "none") +
  theme_default() +
  labs(y = expression(delta^13 * C),
       x = "Date") +
  theme(strip.text = element_text(size = 8, hjust = 0),
        text = element_text(size = 8)) +
  scale_fill_colorblind()

library(patchwork)

isotope_plot = nitrogen_isotope_plot /carbon_isotope_plot

ggsave(isotope_plot, file = "plots/isotope_plot.jpg", width = 6.5, height = 6)
