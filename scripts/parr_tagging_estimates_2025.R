# Purpose: Survival estimates for parr tagging, tag year 2025, migratory year 2026
# Author: Michelle A. Briggs; Ryan N. Kinzer
# Created Date: August 07, 2026


# load libraries
#remotes::install_github("ryankinzer/PITmodelR", ref = "dev", build_vignettes = F, force = T)

library(PITmodelR)
library(tidyverse)
library(marked)

plot_theme <- readRDS("./docs/plot_theme_v2.RDS")

source('./scripts/run_cjs_model_set.R')

# get interrogation site configuration
site_config <- DART_site_config()


##### PARR TAGGING ESTIMATES #######

# get Data ---- 

activity <- 'parr' # used for directory path and hatchery, rst
tag_yr = 2025


# check for mark groups and if interrogation files exist (requires two specific file names "....Mark.rds" and "....INT.rds")

#CTH files are manually downloaded from PTAGIS

file_df <- check_mark_groups_tag_yr(activity = activity, tag_year = 2025)

# we can only run survivals if we have interrogation or complete tag history information
mark_groups <- file_df %>% filter(cth_exists)

no_cth <- file_df %>% filter(!cth_exists)

# Need to map across all the groups contained in mark_groups

mark_df <- map_df(mark_groups$mark_file,
                  ~readRDS(file.path(file_path, .))
)

spp <- c('12W', '11W') #wild spring/summer chinook

#parr tagging locations        
mark_site <- c('BURNLC', 
               'GROUSC', 
               'JOHNSC', 
               'SAEFSF', 
               'SUMITC', 
               'CAPEHC', 
               'CHAMBC', 
               'LAKEC', 
               'MARSHC', 
               'SALRSF', 
               'SECESR',
               'BEAR6C', 
               'HURRIC', 
               'WALLOR', 
               'LOLOC')

# list of tag_ids to use as marks
tag_ids <- mark_df %>%
  filter(species_run_rear_type %in% spp,
         release_site %in% mark_site, #look at what we're filtering here
         !grepl('Y', text_comments), 
         !grepl('Y|M', conditional_comments)) %>% #I think this is removing yearlings and mortalities
  pull(pit_tag)


# Complete Tag History File

cth_df <- map_df(mark_groups$cth_file,
                 ~read_csv(file.path(file_path, .))
) %>%
  select(tag_code = 'Tag Code',
         event_type = 'Event Type Name',
         event_site_name = 'Event Site Name',
         event_datetime = 'Event Date Time Value',
         antenna_id = 'Antenna ID',
         release_site_name = 'Event Release Site Name',
         release_datetime = 'Event Release Date Time Value'
  ) %>%
  filter(tag_code %in% tag_ids) %>%
  mutate(event_site = gsub(" .*$", "", event_site_name),
         release_site = gsub(" .*$", "", release_site_name),
         event_datetime = mdy_hms(event_datetime),
         release_datetime = mdy_hms(release_datetime),
         obs_site = if_else(is.na(release_site),event_site, release_site),
         obs_datetime = if_else(is.na(release_datetime), event_datetime, release_datetime)
  ) %>%
  left_join(site_config, by = c("obs_site", "antenna_id")) %>%
  mutate(date_start = if_else(event_type != 'Observation', obs_datetime, date_start),
         date_end = if_else(event_type != 'Observation', obs_datetime, date_end),
         disposition = if_else(event_type != 'Observation', event_type, disposition),
         censored = if_else(event_type != 'Observation', FALSE, censored)) %>%
  filter(between(obs_datetime, date_start, date_end)) %>%
  arrange(tag_code, obs_datetime) %>%
  group_by(tag_code) %>%
  mutate(
    event = cumsum(obs_site != lag(obs_site, default = first(obs_site))) + 1L
  ) %>%
  ungroup()


ch_df <- cth_df %>%
  group_by(tag_code, obs_site,event_type, event) %>%
  summarise(
    first_obs = min(obs_datetime),
    last_obs  = max(obs_datetime),
    last_disposition = last(disposition[order(obs_datetime)]),
    last_censored = last(censored[order(obs_datetime)]),
    .groups = "drop"
  ) %>%
  arrange(tag_code, first_obs)

# combine mark data with CTH to parse into groups

ch_long <- ch_df %>%
  left_join(mark_df %>%
              select(tag_code = pit_tag,
                     srr = species_run_rear_type,
                     release_site = release_site,
                     mark_release_date = release_date,
                     release_season = event_season,
                     brood_year,
                     migration_year,
                     text_comments,
                     conditional_comments),
            by = 'tag_code')
# 
# %>%
#   filter(release_site == 'JOHTRP',
#          release_season == 'Summer/Fall')


#check some numbers
table(ch_long$last_disposition)
table(ch_long$obs_site)
table(ch_long$event_type)
table(ch_long$release_site)
table(ch_long$release_season)
table(ch_long$srr)
table(ch_long$text_comments) 
table(ch_long$conditional_comments) 


lgr <- c('GRJ', 'GRS') #s = spillway, #j = juvenile bypass
down <- c("GOJ", "LMJ", "ICH", "MCJ", "JDJ", "B2J", "BCC",
          "PD5", "PD6", "PD7", "PD8", "PDO", "PDW", "TWX")

obs_tbl <- tibble(
  release_site = c('BURNLC', 
                   'GROUSC', 
                   'JOHNSC', 
                   'SAEFSF', 
                   'SUMITC', 
                   'CAPEHC', 
                   'CHAMBC', 
                   'LAKEC', 
                   'MARSHC', 
                   'SALRSF', 
                   'SECESR',
                   'BEAR6C', 
                   'HURRIC', 
                   'WALLOR', 
                   'LOLOC'),
  obs_loc = list(
    list(
      BURNLC = "BURNLC",
      ESS = "ESS",
      SFG = "SFG",
      LGR = lgr,
      Down = down
    ), 
    list(
      GROUSC = "GROUSC",
      ZEN = "ZEN", 
      SFG = "SFG", 
      LGR = lgr, 
      Down = down
    ), 
    list(
      JOHNSC = "JOHNSC",
      ESS = "ESS", 
      SFG = "SFG", 
      LGR = lgr, 
      Down = down
    ), 
    list(
      SAEFSF = "SAEFSF",
      ESS = "ESS", 
      SFG = "SFG", 
      LGR = lgr, 
      Down = down
    ), 
    list(
      SUMITC = "SUMITC",
      ZEN = "ZEN", 
      SFG = "SFG", 
      LGR = lgr, 
      Down = down
    ), 
    list(
      CAPEHC = "CAPEHC",
      MAR = "MAR", 
      LGR = lgr, 
      Down = down
    ), 
    list(
      CHAMBC = "CHAMBC",
      LGR = lgr, 
      Down = down
    ), 
    list(
      LAKEC = "LAKEC",
      ZEN = "ZEN", 
      SFG = "SFG", 
      LGR = lgr, 
      Down = down
    ), 
    list(
      MARSHC = "MARSHC",
      MAR = "MAR", 
      LGR = lgr, 
      Down = down
    ), 
    list(
      SALRSF = "SALRSF",
      KRS = "KRS", 
      SFG = "SFG", 
      LGR = lgr, 
      Down = down
    ), 
    list(
      SECESR = "SECESR",
      ZEN = "ZEN", 
      SFG = "SFG", 
      LGR = lgr, 
      Down = down
    ), 
    list(
      BEAR6C = "BEAR6C",
      WR2 = "WR2", 
      WR1 = "WR1", 
      LGR = lgr, 
      Down = down
    ),
    list(
      HURRIC = "HURRIC", 
      WR2 = "WR2", 
      WR1 = "WR1", 
      LGR = lgr, 
      Down = down
    ),  
    list(
      WALLOR = "WALLOR", 
      WR2 = "WR2", 
      WR1 = "WR1", 
      LGR = lgr, 
      Down = down
    ), 
    list(
      LOLOC = "LOLOC", 
      LC2 = "LC2",
      LC1 = "LC1", 
      CWR = "CWR", 
      LGR = lgr, 
      Down = down
    )
  )
)

ch_long_list <- ch_long %>%
  group_by(srr, release_site, release_season) %>%
  nest() %>%
  left_join(obs_tbl)

ch_long %>%
  group_by(srr, release_site, release_season) %>% 
  count()  

parr_estimates <- run_cjs_model_set(ch_long_list)

parr_survival_estimates <- parr_estimates$estimates
parr_survival_tbl <- parr_estimates$model_tbl


saveRDS(ch_long, file = "data/parr_tagging/ch_parr_2025.RDS")
saveRDS(parr_estimates, file = "data/parr_tagging/parr_estimates_2025.RDS")
####FIGURES#####


parr_survival_estimates %>%
  filter(metric == "survival", interval != 4) %>%
  unite(col = "int_label", from_site, to_site, sep = " to ", remove = FALSE) %>%
  ggplot(aes(x = as.factor(interval), y = cjs_est)) + 
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = marray_lcl, ymax = marray_ucl), width = 0.1) +
  plot_theme + 
  facet_wrap(~release_site) +
  scale_y_continuous(limits = c(0, 1)) + 
  theme(strip.background = element_rect(fill = "grey", color = "black")) + 
  ylab("Survival") + 
  xlab("Interval")

parr_survival_estimates %>%
  filter(metric == "detection", interval %in% c(2, 3, 4)) %>%
  ggplot(aes(x = as.factor(interval), y = marray_est)) + 
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = marray_lcl, ymax = marray_ucl), width = 0.1) +
  plot_theme + 
  facet_wrap(~release_site) + 
  scale_y_continuous(limits = c(0, 1)) + 
  theme(strip.background = element_rect(fill = "grey", color = "black")) + 
  ylab("Detection probability") + 
  xlab("Site")

parr_estimates_long <- sf_parr_survival_estimates %>%
  pivot_longer(cols = matches("^(cjs|mscjs|marray)_(est|lcl|ucl)$"),
               names_to = c("model", ".value"),
               names_pattern = "(cjs|mscjs|marray)_(est|lcl|ucl)")

parr_estimates_long %>%
  filter(metric == "survival", interval != 4) %>%
  ggplot(aes(x = as.factor(interval), y = est, color = model)) + 
  geom_point(size = 3, position = position_dodge(width = 0.25)) + 
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.1, linewidth = 1,
                position = position_dodge(width = 0.25)) +
  plot_theme + 
  facet_wrap(~release_site) + 
  scale_y_continuous(limits = c(0, 1)) + 
  theme(strip.background = element_rect(fill = "white", color = "black")) + 
  scale_color_brewer(palette = "Dark2") + 
  ylab("Survival") + 
  xlab("Interval")
