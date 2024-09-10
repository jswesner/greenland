library(tidyverse)
library(tidybayes)
library(ggthemes)
library(brms)

# load data
thule_birds_clean = read_csv(file = "data/thule_birds_clean.csv") %>% 
  mutate(mean_julian = mean(julian),
         sd_julian = sd(julian),
         julian_z = (julian - mean_julian)/sd_julian,
         hg_mean = hg/mean(hg)) %>% 
  group_by(species_common) %>% 
  mutate(mean_hg = mean(hg),
         hg_mean_species = hg/mean_hg)

# load model
brm_smooth_n_mi = readRDS(file = "models/brm_smooth_n_mi.rds")
brm_smooth_n = readRDS(file = "models/brm_smooth_n.rds")

# make data grid to predict over
# Note: Atlantic Puffin don't have julian_z > ~ 0.75, so their recent estimates are model predictions

species_means = thule_birds_clean %>% distinct(species_common, mean_hg)
mean_julian = mean(thule_birds_clean$julian)
sd_julian = sd(thule_birds_clean$julian)


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
  group_seqs[[i]] = tibble(julian_z = seq(min(temp[[i]]$julian_z),
                                          max(temp[[i]]$julian_z),
                                          length.out = 30)) %>% 
    mutate(species_common = unique(temp[[i]]$species_common))
}


post_preds = bind_rows(group_seqs) %>%
  mutate(nitrogen_s = 0) %>% 
  add_epred_draws(brm_smooth_n, re_formula = NULL) %>% 
  left_join(species_means) %>% 
  mutate(julian = (julian_z*sd_julian) + mean_julian,
         date = as.Date(julian)) %>% 
  group_by(species_common) %>% 
  mutate(median = median(.epred*mean_hg))


post_timeseries = post_preds %>% 
  ggplot(aes(x = julian_z, y = .epred*mean_hg, fill = species_common)) + 
  stat_lineribbon(alpha = 0.4, .width = 0.95, linewidth = 0.1) + 
  facet_wrap(~reorder(species_common, -median), nrow = 1) +
  geom_point(data = thule_birds_clean %>% left_join(post_preds %>% ungroup %>% 
                                                      distinct(species_common, median)), 
             aes(y = hg, color = species_common),
             shape = 1,
             size = 0.1) +
  guides(fill = "none",
         color = "none") +
  theme_default() +
  labs(y = "Hg (ng/gWM blood)",
       x = "Date") +
  theme(strip.text = element_text(size = 11),
        text = element_text(size = 11)) +
  scale_fill_colorblind() + 
  scale_color_colorblind()

ggsave(post_timeseries, file = "plots/post_timeseries.jpg", width = 6.5, height = 3, dpi = 400)

post_means = post_preds %>% 
  group_by(species_common, .draw, median) %>% 
  mutate(hg = .epred*mean_hg) %>% 
  reframe(hg = mean(hg))

post_mean_plot = post_means %>% 
  ggplot(aes(x = reorder(species_common, median),
             y = hg)) +
  geom_jitter(data = thule_birds_clean %>% left_join(post_preds %>% ungroup %>% 
                                                       distinct(species_common, median)), 
              aes(y = hg),
              width = 0.2,
              size = 0.1,
              shape = 1) + 
  geom_boxplot(aes(group = species_common, fill = species_common),
               outlier.shape = NA,
               linewidth = 0.1) +
  guides(fill = "none",
         color = "none") +
  theme_default() +
  labs(y = "Hg (ng/g blood)",
       x = "Species") +
  theme(strip.text = element_text(size = 11),
        text = element_text(size = 11)) +
  scale_fill_colorblind() + 
  scale_color_colorblind() +
  geom_hline(yintercept = c(90, 600, 1300))

ggsave(post_mean_plot, file = "plots/post_mean_plot.jpg", width = 5, height = 5, dpi = 400)


# plot nitrogen 
post_n = tibble(nitrogen_s = seq(min(brm_smooth_n$data$nitrogen_s),
                        max(brm_smooth_n$data$nitrogen_s),
                        length.out = 30)) %>% 
  expand_grid(species_common = unique(brm_smooth_n$data$species_common)) %>% 
  mutate(julian_z = 0) %>% 
  add_epred_draws(brm_smooth_n, file = "models/brm_smooth_n.rds", re_formula = NULL) %>% 
  left_join(species_means)


post_n %>% 
  ggplot(aes(x = nitrogen_s, y = .epred*mean_hg)) + 
  stat_lineribbon(.width = 0.95) +
  geom_point(data = thule_birds_clean, aes(y = hg)) +
  facet_wrap(~species_common) +
  NULL
