library(tidyverse)
library(tidybayes)
library(ggthemes)
library(brms)
theme_set(brms::theme_default())

# load data
thule_birds_clean = read_csv(file = "data/thule_birds_clean.csv") %>% 
  mutate(species_abbreviation = case_when(species_common == "Dovekie" ~ "e) DOVE",
                                          species_common == "Black-legged Kittiwake" ~ "d) BLKI",
                                          species_common == "Atlantic Puffin" ~ "c) ATPU",
                                          species_common == "Black Guillemot" ~ "b) BLGU",
                                          species_common == "Thick-billed Murre" ~ "a) TBMU"))

# load model
brm_smooth_nc_gamma = readRDS("models/brm_smooth_nc_gamma.rds")

# load posteriors
post_preds_hg = readRDS(file = "posteriors/post_preds_hg.rds") 

post_means = post_preds_hg %>% 
  group_by(species_common, .draw) %>% 
  reframe(hg = mean(hg))

post_c_overall = readRDS(file = "posteriors/post_c_overall.rds")


# summarize ---------------------------------------------------------------
# raw summary
raw_table = thule_birds_clean %>% 
  group_by(species_common) %>% 
  add_tally() %>% 
  reframe(min = min(hg),
          max = max(hg),
          n = n,
          mean_n = mean(nitrogen, na.rm = T),
          mean_c = mean(carbon, na.rm = T)) %>% 
  distinct() %>% 
  mutate_if(is.numeric, round, 1) %>% 
  mutate(min = round(min, 0),
         max = round(max, 0)) %>% 
  mutate(hg_range = paste0(min, " to ", max)) %>% 
  select(-min, -max)

# species median Hg
post_preds_hg %>% 
  group_by(species_common) %>% 
  reframe(mean = mean(hg),
          sd = sd(hg)) %>% 
  arrange(mean)

table_2 = post_preds_hg %>% 
  group_by(species_common) %>% 
  reframe(mean = mean(hg),
          sd = sd(hg)) %>% 
  mutate(mean = round(mean, 0),
         sd = round(sd, 0)) %>% 
  arrange(mean) %>% 
  mutate(mean_sd = paste(mean, "pm", sd)) %>% 
  left_join(raw_table) %>% 
  select(-mean, -sd)

write_csv(table_2, file = "tables/table_2.csv")

# probability of positive slopes between Hg and N
tibble(nitrogen_s = c(0, 1)) %>% 
  expand_grid(species_common = unique(thule_birds_clean$species_common)) %>% 
  expand_grid(carbon_s = c(-1, 0, 1)) %>% 
  expand_grid(julian_z = c(-2, -1, 0, 1, 2)) %>% 
  add_epred_draws(brm_smooth_nc_gamma) %>%
  left_join(species_means) %>% 
  mutate(.epred = .epred*mean_hg) %>% 
  # group_by(nitrogen_s, .draw) %>% 
  group_by(nitrogen_s, .draw, species_common) %>% 
  reframe(.epred = mean(.epred)) %>% 
  pivot_wider(names_from = nitrogen_s, values_from = .epred) %>% 
  mutate(slope = (`1` - `0`)/`0`) %>% 
  group_by(species_common) %>%
  reframe(prob_pos = sum(slope>0)/max(.draw))
  median_qi(slope)

# probability of positive slopes between Hg and C
tibble(carbon_s = c(0, 1)) %>% 
  expand_grid(species_common = unique(thule_birds_clean$species_common)) %>% 
  expand_grid(nitrogen_s = c(-1, 0, 1)) %>% 
  expand_grid(julian_z = c(-2, -1, 0, 1, 2)) %>% 
  add_epred_draws(brm_smooth_nc_gamma) %>%
  left_join(species_means) %>% 
  mutate(.epred = .epred*mean_hg) %>% 
  group_by(carbon_s, .draw) %>% 
  reframe(.epred = mean(.epred)) %>% 
  pivot_wider(names_from = carbon_s, values_from = .epred) %>% 
  mutate(slope = (`1` - `0`)/`0`) %>% 
  # group_by(species_common) %>%
  reframe(prob_pos = sum(slope>0)/max(.draw))
  median_qi(slope)


# dates with min and max hg
post_preds_hg %>% 
  group_by(julian_z, date, species_common) %>%
  median_qi(hg) %>% 
  group_by(species_common) %>% 
  filter(hg == max(hg)|hg == min(hg)) %>% 
  arrange(species_common, hg)

# difference between first and last year mean Hg
post_last_minus_first = post_preds_hg %>% 
  group_by(species_common) %>% 
  mutate(year = year(date)) %>% 
  filter(year == max(year)|year == min(year)) %>% 
  group_by(year, species_common, .draw) %>% 
  reframe(mean_hg = mean(hg)) %>%
  group_by(species_common) %>% 
  mutate(first_last = case_when(year == min(year) ~ "first",
                                year == max(year) ~ "last")) %>% 
  select(-year) %>% 
  pivot_wider(names_from = first_last, values_from = mean_hg) %>% 
  mutate(diff_year = (last/first) - 1) 
  
post_last_minus_first %>% 
  group_by(species_common) %>%
  pivot_longer(cols = c(first, last)) %>% 
  group_by(name, species_common) %>% 
  median_qi(value)

post_last_minus_first %>% 
  group_by(species_common) %>% 
  median_qi(diff_year)

post_last_minus_first %>% 
  group_by(species_common) %>% 
  reframe(prob_higher = sum(diff_year>1)/max(.draw))



# isotope timeseries ------------------------------------------------------
isotope_data = readRDS(file = "data/isotope_data.rds")
brm_smooth_isotope_model = readRDS(file = "models/brm_smooth_isotope_model.rds")

linear_julian_seq = seq(min(isotope_data$julian_z), 
                        max(isotope_data$julian_z),
                        length.out = 100)

post_preds_isotopes = tibble(julian_z = round(quantile(linear_julian_seq), 2)) %>% 
  expand_grid(species_isotope = unique(brm_smooth_isotope_model$data$species_isotope)) %>% 
  left_join(isotope_data %>% ungroup %>% distinct(species_common, isotope, species_isotope, mean_julian, sd_julian, isotope_mean)) %>% 
  mutate(date = as.Date(julian_z*sd_julian) + mean_julian,
         year = year(date)) %>% 
  add_epred_draws(brm_smooth_isotope_model) 
  
isotope_posts_by_year = post_preds_isotopes %>% 
  group_by(species_common, isotope, year, .draw) %>% 
  reframe(.epred = mean(.epred*isotope_mean)) %>% 
  pivot_wider(names_from = year, values_from = .epred) %>% 
  mutate(first_last = `2023` - `2010`,
         first_middle = `2017` - `2010`,
         middle_last = `2023` - `2017`) %>% 
  pivot_longer(cols = c(first_last, first_middle, middle_last))

isotope_posts_by_year %>% 
  group_by(species_common, isotope, name) %>% 
  filter(name == "first_last") %>% 
  # reframe(prob_lower = sum(value<0)/max(.draw))
  median_qi(value) 
  
isotope_posts_by_year %>% 
  group_by(isotope, name) %>% 
  filter(name == "first_last") %>% 
  group_by(species_common, isotope, .draw) %>% 
  reframe(value = mean(value)) %>% 
  # group_by(isotope, .draw) %>% 
  # reframe(value = mean(value)) %>% 
  group_by(isotope, species_common) %>%
  # reframe(prob_diff = sum(value>0)/max(.draw))
  median_qi(value) %>% 
  ggplot(aes(x = value, xmin = .lower, xmax = .upper, y = species_common, color = isotope)) +
  geom_pointrange()

  
readRDS("posteriors/post_isotopes.rds") %>% 
  # filter(species_common == "Dovekie") %>% 
  group_by(species_common, isotope) %>% 
  median_qi(.epred)


# hg and isotopes together -------------------------------------------------

# last minus first table
post_isotope_simple = post_preds_isotopes %>% ungroup %>% 
  left_join(isotope_data %>% ungroup %>% distinct(species_common, isotope, species_isotope, mean_julian, sd_julian, isotope_mean)) %>%
  mutate(.epred = .epred*isotope_mean) %>% 
  select(date, species_common, isotope, .epred, .draw) %>% 
  rename(response = isotope) 

post_hg_simple = post_preds_hg %>% select(date, species_common, hg, .draw) %>% rename(.epred = hg) %>% 
  mutate(response = "hg")

post_all = post_isotope_simple %>% bind_rows(post_hg_simple) %>% mutate(year = year(date))

post_all_years = post_all %>% 
  group_by(species_common, response, year, .draw) %>% 
  reframe(.epred = mean(.epred)) 
  
last_minus_first_table = post_all_years %>% 
  group_by(species_common, response) %>% 
  filter(year == min(year)|year == max(year)) %>% 
  pivot_wider(names_from = year, values_from = .epred) %>% 
  mutate(diff = case_when(is.na(`2023`) ~ `2019` - `2010`,
                          TRUE ~ `2023` - `2010`)) %>% 
  group_by(species_common, response) %>% 
  reframe(median = median(diff),
          .lower = quantile(diff, probs = 0.025),
          .upper = quantile(diff, probs = 0.975),
          prob_higher = sum(diff < 0)/max(.draw)) %>% 
  mutate(response = as.factor(response),
         response = fct_relevel(response, "hg", "nitrogen")) %>% 
  arrange(response) %>% 
  mutate_if(is.numeric, round, 1) %>% 
  mutate(comparison = "last - first") %>% 
  mutate(slope_summary = paste0(median, " (", .lower, "±", .upper, ")")) %>% 
  select(species_common, response, comparison, slope_summary)

write_csv(last_minus_first_table, file = "tables/last_minus_first_table.csv")


# slopes table
# probability of positive slopes between Hg and N
hg_n_slopes = tibble(nitrogen_s = c(0, 1)) %>% 
  expand_grid(species_common = unique(thule_birds_clean$species_common)) %>% 
  expand_grid(carbon_s = c(-1, 0, 1)) %>% 
  expand_grid(julian_z = c(-2, -1, 0, 1, 2)) %>% 
  add_epred_draws(brm_smooth_nc_gamma) %>%
  left_join(species_means) %>% 
  mutate(.epred = .epred*mean_hg) %>% 
  # group_by(nitrogen_s, .draw) %>% 
  group_by(nitrogen_s, .draw, species_common) %>% 
  reframe(.epred = mean(.epred)) %>% 
  pivot_wider(names_from = nitrogen_s, values_from = .epred) %>% 
  mutate(slope = (`1` - `0`)/`0`) %>% 
  group_by(species_common) %>%
  median_qi(slope) %>% 
  mutate(response = "hg",
         predictor = "deltaN")

# probability of positive slopes between Hg and C
hg_c_slopes = tibble(carbon_s = c(0, 1)) %>% 
  expand_grid(species_common = unique(thule_birds_clean$species_common)) %>% 
  expand_grid(nitrogen_s = c(-1, 0, 1)) %>% 
  expand_grid(julian_z = c(-2, -1, 0, 1, 2)) %>% 
  add_epred_draws(brm_smooth_nc_gamma) %>%
  left_join(species_means) %>% 
  mutate(.epred = .epred*mean_hg) %>% 
  group_by(carbon_s, .draw, species_common) %>% 
  reframe(.epred = mean(.epred)) %>% 
  pivot_wider(names_from = carbon_s, values_from = .epred) %>% 
  mutate(slope = (`1` - `0`)/`0`) %>% 
  group_by(species_common) %>%
  median_qi(slope) %>% 
  mutate(response = "hg",
         predictor = "deltaC")

slope_hg_n_c_table = bind_rows(hg_n_slopes, hg_c_slopes) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  mutate(slope_summary = paste0(slope, " (", .lower, "±", .upper, ")")) %>% 
  select(species_common, response, predictor, slope_summary) 

write_csv(slope_hg_n_c_table, file = "tables/slope_hg_n_c_table.csv")
