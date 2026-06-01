# Purpose: Simulate capture history data and evaluate estimates.
# Author: Ryan N. Kinzer
# Created Date: April 28, 2026


# load libraries

dnload_new_version = TRUE

if (dnload_new_version) {
  
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  
  remotes::install_github(
    "ryankinzer/PITmodelR",
    ref = "dev"
  )
}

library(PITmodelR)
library(tidyverse)

source('./scripts/run_cjs_model_set.R')


# sim function
run_one_sim <- function(iter,
                        n_init = 1000,
                        n_add = c(0, 0, 0),
                        S = c(0.5, 0.6),
                        p = c(1, 0.2, 0.05),
                        remove_prob = c(0.0, 0.0, 0.0)) {
  
  sites <- paste0("R", 1:length(p))
  
  
  sim_events <- simulate_cjs_events(
    n_init = n_init,
    n_add = n_add,
    sites = sites,
    S = S,
    p = p,
    remove_prob = remove_prob,
    seed = iter
  )
  
  obs_loc <- as.list(sites)
  names(obs_loc) <- sites
  
  sim_list <- tibble(
    sim = iter,
    srr = "sim",
    release_site = sites[1],
    release_season = "test",
    data = list(sim_events),
    obs_loc = list(obs_loc)
  )
  
  fit <- run_cjs_model_set(sim_list)
  
  fit$estimates %>%
    mutate(sim = iter, .before = 1)
}


# generated simulation results

survival <- c(0.5, 0.6)
detection <- c(1, 0.2, 0.02)
removal <- c(0.0, 0.0, 0.0)


sim_results <- map_dfr(
  1:100,
  ~ run_one_sim(.x,
                S = survival,
                p = detection,
                remove_prob = removal)
)

true_vals <- tibble(
  metric = c("survival", "survival", "detection", "detection"),
  interval = c(1, 2, 2, 3),
  parameter = c("S_1", "S_2", "p_2", "p_3"),
  true = c(0.5, 0.6, 0.2, 0.05)
)

sim_ests <- sim_results %>%
  left_join(true_vals, by = c("metric", "interval")) %>%
  select(
    sim, metric, interval, from_site, to_site, true,
    cjs_est, cjs_lcl, cjs_ucl,
    mscjs_est, mscjs_lcl, mscjs_ucl,
    marray_est, marray_lcl, marray_ucl
  ) %>%
  pivot_longer(
    cols = matches("^(cjs|mscjs|marray)_(est|lcl|ucl)$"),
    names_to = c("model", ".value"),
    names_pattern = "(cjs|mscjs|marray)_(est|lcl|ucl)"
  ) %>%
  filter(!is.na(true), !is.na(est))


sim_summary <- sim_ests %>%
  group_by(model, metric, interval, from_site, to_site) %>%
  summarise(
    n = n(),
    true = first(true),
    mean_est = mean(est, na.rm = TRUE),
    bias = mean(est - true, na.rm = TRUE),
    percent_bias = 100 * mean((est - true) / true, na.rm = TRUE),
    mse = mean((est - true)^2, na.rm = TRUE),
    rmse = sqrt(mse),
    coverage = mean(lcl <= true & ucl >= true, na.rm = TRUE),
    mean_ci_width = mean(ucl - lcl, na.rm = TRUE),
    .groups = "drop"
  )

sim_summary

sim_ests %>%
  filter(interval == 1,
         metric == 'survival') %>%
  ggplot(aes(x = est, fill = model)) +
  geom_histogram() +
  facet_wrap(~model)


#######



obs_loc = list(
  R1 = "R1",
  R2 = "R2",
  R3 = "R3")

sim_events <- simulate_cjs_events(
  n_init = 1000,   # number of marks originally released
  n_add = c(0, 0, 0), # marks added at each detection/release site
  sites = paste0('R', 1:length(detection)), # release sites
  S = survival, # true survivals probs. (1 less than number of sites with the first value being the survival from R1 to R2)
  p = detection, # true detection probs.
  remove_prob = removal, #randomly select marks that are detected and then removed at each site
  seed = 123
)

sim_list <- tibble(
  srr = "sim",
  release_site = "R1",
  release_season = "test",
  data = list(sim_events),
  obs_loc = list(obs_loc)
)

tmp <- run_cjs_model_set(sim_list)

true_vals <- tibble(
  metric = c("survival", "survival", "detection", "detection", "detection"),
  interval = c(1, 2, 1, 2, 3),
  parameter = c("S_1", "S_2", "p_1", "p_2", "p_3"),
  true = c(survival, detection)
)

ests <- tmp$estimates %>%
  left_join(true_vals, by = c("metric", "interval")) %>%
  filter(
    interval == 1,
    metric == "survival"
  ) %>%
  select(true, contains("_est")) %>%
  pivot_longer(
    cols = contains("_est"),
    names_to = "model",
    values_to = "est"
  ) %>%
  mutate(
    model = gsub("_est", "", model),
    bias = true,
    percent_bias = 100 * (est-true) / true
  )




ch_list <- build_capture_histories(
  tag_history = sim_events,
  locs_def = obs_locs,
  site_col = "obs_site",
  tag_col = "tag_code",
  time_col = "first_obs",
  censor_col = "last_censored"
)

ch_data <- ch_list$ch_data
ch_data$ch <- gsub("2", "1", ch_data$ch)
m_array <- ch_list$m_array
ms_data <- build_multistate_histories(ch_data)

fit_cjs <- fit_marked_cjs(ch_data)

fit_mscjs <- fit_marked_mscjs(ms_data)

fit_marray <- fit_marray_cjs(m_array)

fit_cjs$phi
fit_mscjs$phi
fit_marray$phi

fit_cjs$p
fit_mscjs$p
fit_marray$p

fit_marray$lambda

true_values <- data.frame(
  parameter = c(
    paste0("S_", 1:4),
    paste0("p_", 2:5)
  ),
  true = c(
    c(0.85, 0.75, 0.90, 0.80),
    c(0.40, 0.50, 0.60, 0.70)
  )
)

true_values