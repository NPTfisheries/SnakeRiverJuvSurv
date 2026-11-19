# Purpose: Survival estimates for RST tagging for comparison with parr tagging estimates
# 29 October 2026

#Estimates for JOHTRP 2023 - 2025

# load libraries
#remotes::install_github("ryankinzer/PITmodelR", ref = "dev", build_vignettes = F, force = T)

library(PITmodelR)
library(tidyverse)
library(marked)
library(RTMB)

plot_theme <- readRDS("./docs/plot_theme_v2.RDS")

source('./scripts/run_cjs_model_set.R')

# get interrogation site configuration
site_config <- DART_site_config()


##### RST TAGGING ESTIMATES #######

#produce estimates for Johnson Creek RST fish tagged in July and August, 2021 - 2024

# get data

activity <- 'rst' # used for directory path

year <- seq(2021, 2025, by = 1)


#Define observation sites
#same for all years

lgr <- c('GRJ', 'GRS', 'LGRRRR') #s = spillway, #j = juvenile bypass
down <- c("GOJ", "LMJ", "ICH", "MCJ", "JDJ", "B2J", "BCC",
          "PD5", "PD6", "PD7", "PD8", "PDO", "PDW", "TWX")

obs_tbl <- tibble(
  release_site = c('JOHTRP'),
  obs_loc = list(
    list(
      JOHTRP = "JOHTRP",
      ESS = "ESS", 
      SFG = "SFG", 
      LGR = lgr, 
      Down = down
    )
  )
)

#Johnson Creek didn't run past June 30 in 2025

for (y in 1:length(year)) {

tag_yr = year[y]


# check for mark groups and if interrogation files exist (requires two specific file names "....Mark.rds" and "....INT.rds")

#CTH files are manually downloaded from PTAGIS

file_df <- check_mark_groups_tag_yr(activity = activity, tag_year = tag_yr)

# we can only run survivals if we have interrogation or complete tag history information
mark_groups <- file_df %>% filter(cth_exists)

no_cth <- file_df %>% filter(!cth_exists)

# Need to map across all the groups contained in mark_groups


mark_df <- map_df(mark_groups$mark_file,
                  ~readRDS(file.path(file_path, .))
)

#filter for summer tagging in July and August

ggplot(mark_df, aes(x = event_date)) +
  geom_histogram() +
  facet_wrap(~release_site, scale = "free_y")

mark_df <- mark_df %>%
  filter(event_date <= paste0(tag_yr, '-08-31') & event_date >= paste0(tag_yr, '-07-01'))

mark_by_site <- mark_df %>% 
  group_by(release_site) %>%
  count()

spp <- c('12W', '11W') #wild spring/summer chinook

#parr tagging locations        
mark_site <- c('JOHTRP')

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
                     #release_season = event_season,
                     #mark_season,
                     brood_year,
                     migration_year,
                     text_comments,
                     conditional_comments),
            by = 'tag_code')
  
#save year-specific ch for checking
 
assign(paste0("ch_long_", tag_yr), ch_long)


#check some numbers
# table(ch_long$last_disposition)
# table(ch_long$obs_site)
# table(ch_long$event_type)
# table(ch_long$release_site)
# table(ch_long$release_season)
# table(ch_long$srr)
# table(ch_long$mark_season)
# table(ch_long$text_comments) 
# table(ch_long$conditional_comments) 


ch_long_list <- ch_long %>%
  group_by(srr, release_site) %>%
  nest() %>%
  left_join(obs_tbl)



#define all as summer chinook so they will run together
ch_long_list <- ch_long_list %>%
  mutate(srr = '12W')

ch_long_list <- ch_long_list %>%
  mutate(
    ch = map2(
      data,
      obs_loc,
      ~ build_capture_histories(
        tag_history = .x,
        locs_def = .y,
        site_col = "obs_site",
        tag_col = "tag_code",
        time_col = "first_obs",
        censor_col = "last_censored",
        covariate_cols = NULL
      )
    )
  )

#individual m array for each release group
marray_list <- ch_long_list %>%
  mutate(m_array = map(
    ch, ~(marray = .x$m_array))) %>%
  select(srr, release_site, m_array) %>%
  rename(release_group = release_site)

#full shared detection site model
shared_p <- "All"
shared_phi <- "None"

fit_multiple <- fit_marray_cjs_multiple(marray_list, 
                                              shared_p, shared_phi)


saveRDS(fit_multiple, file = paste0("output/estimates/summer_JOHTRP_", tag_yr, ".RDS"))

}


################ PARR TAGGING ESTIMATES ######################

#from NOAA and ODFW
#project codes OGR and GAA

#check on catherine creek, grande ronde, and imnaha IPTDS sites

lgr <- c('GRJ', 'GRS', 'LGRRRR') #s = spillway, #j = juvenile bypass
down <- c("GOJ", "LMJ", "ICH", "MCJ", "JDJ", "B2J", "BCC",
          "PD5", "PD6", "PD7", "PD8", "PDO", "PDW", "TWX")
SC12 <- c('SC1', 'SC2') #pool detections at SC1 and SC2
LC12 <- c('LC1', 'LC2') #pool detections at LC1 and LC2
IR12 <- c('IR1', 'IR2')


obs_tbl <- tibble(
  release_site = c('AMERR', 
                   'BEAR6C',
                   'BEARVC',
                   'BIG2C',
                   'BURNLC',
                   'CAPEHC', 
                   'CATCNF', 
                   'CATCSF',
                   'CATHEC',
                   'CHAMBC', 
                   'CROOKR', 
                   'ELKC',
                   'GRAND2',
                   'GROUSC', 
                   'HURRIC', 
                   'IMNAHR',
                   'JOHNSC', 
                   'KNOXB',
                   'LAKEC',  
                   'LOLOC',
                   'LOONC',
                   'LOSTIR',
                   'MARSHC', 
                   'MINAMR',
                   'NEWSOC', 
                   'REDR',
                   'SAEFSF', 
                   'SALRSF', 
                   'SECESR',
                   'SUMITC', 
                   'VALEYC',
                   'WALLOR' 
  ),
  obs_loc = list(
    list(
      AMERR = "AMERR", 
      SC4 = "SC4",
      SC3 = "SC3", 
      SC12 = SC12, 
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
      BEARVC = "BEARVC", 
      BRC = "BRC",
      LGR = lgr,
      Down = down
    ),
    list(
      BIG2C = "BIG2C", 
      TAY = "TAY", 
      LGR = lgr,
      Down = down
    ),
    list(
      BURNLC = "BURNLC",
      ESS = "ESS",
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
      CATCNF = "CATCNF",
      UGR = "UGR", 
      LGR = lgr, 
      Down = down
    ),
    list(
      CATCSF = "CATCSF",
      UGR = "UGR", 
      LGR = lgr, 
      Down = down
    ),
    list(
      CATHEC = "CATHEC",
      UGR = "UGR", 
      LGR = lgr, 
      Down = down
    ),
    list(
      CHAMBC = "CHAMBC",
      LGR = lgr, 
      Down = down
    ), 
    list(
      CROOKR = "CROOKR",
      CRA = "CRA",
      SC4 = "SC4",
      SC3 = "SC3", 
      SC12 = SC12,
      LGR = lgr, 
      Down = down
    ), 
    list(
      ELKC = "ELKC",
      LGR = lgr, 
      Down = down
    ),
    list(
      GRAND2 = "GRAND2",
      UGS = "UGS", 
      UGR = "UGR", 
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
      HURRIC = "HURRIC", 
      WR2 = "WR2", 
      WR1 = "WR1", 
      LGR = lgr, 
      Down = down
    ),  
    list(
      IMNAHR = "IMNAHR",
      IR3 = "IR3",
      IR12 = IR12, 
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
      KNOXB = "KNOXB", 
      KRS = "KRS", 
      SFG = "SFG", 
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
      LOLOC = "LOLOC", 
      LC12 = LC12,
      LGR = lgr, 
      Down = down
    ), 
    list(
      LOONC = "LOONC", 
      LGR = lgr, 
      Down = down
    ),
    list(
      LOSTIR = "LOSTIR",
      WR2 = "WR2",
      WR1 = "WR1",
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
      MINAMR = "MINAMR",
      MR1 = "MR1",
      WR1 = "WR1",
      LGR = lgr,
      Down = down
    ),
    list(
      NEWSOC = "NEWSOC",
      SC4 = "SC4",
      SC3 = "SC3", 
      SC12 = SC12,
      LGR = lgr, 
      Down = down
    ), 
    list(
      REDR = "REDR",
      SC4 = "SC4",
      SC3 = "SC3", 
      SC12 = SC12,
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
      SUMITC = "SUMITC",
      ZEN = "ZEN", 
      SFG = "SFG", 
      LGR = lgr, 
      Down = down
    ), 
    list(
      VALEYC = "VALEYC",
      VC2 = "VC2",
      LGR = lgr,
      Down = down
    ),
    list(
      WALLOR = "WALLOR", 
      WR2 = "WR2", 
      WR1 = "WR1", 
      LGR = lgr, 
      Down = down
    )
  )
)


activity <- "parr"

tag_yr <- 2025

#run this for each year separately

# check for mark groups and if interrogation files exist (requires two specific file names "....Mark.rds" and "....INT.rds")

#CTH files are manually downloaded from PTAGIS

file_df <- check_mark_groups_tag_yr(activity = activity, tag_year = tag_yr)

# we can only run survivals if we have interrogation or complete tag history information
mark_groups <- file_df %>% filter(cth_exists)

no_cth <- file_df %>% filter(!cth_exists)

# Need to map across all the groups contained in mark_groups

mark_df <- map_df(mark_groups$mark_file,
                  ~readRDS(file.path(file_path, .))
)

#remove groups with fewer than 50 tags (more on the conservative side)

ggplot(mark_df, aes(x = event_date)) +
  geom_histogram() +
  facet_wrap(~release_site, scale = "free_y")

#filter for summer tagging only
mark_df <- mark_df %>%
  filter(event_date <= paste0(tag_yr, '-08-31') & event_date >= paste0(tag_yr, '-07-01'))

mark_by_site <- mark_df %>% 
  group_by(release_site) %>%
  count() %>%
  filter(n >= 50)

sample_size <- mark_by_site$release_site

spp <- c('12W', '11W') #wild spring/summer chinook

#could specify this manually, but for now we'll use what's in the dataset     
mark_site <- unique(mark_df$release_site)


# list of tag_ids to use as marks
#filter by species, adequate sample site, in listed mark site, removing yearling and morts
tag_ids <- mark_df %>%
  filter(species_run_rear_type %in% spp,
         release_site %in% sample_size, 
         !grepl('Y', text_comments), 
         !grepl('Y|M', conditional_comments)) %>% #I think this is removing yearlings and mortalities
  pull(pit_tag)


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
                     #release_season = event_season,
                     #mark_season,
                     brood_year,
                     migration_year,
                     text_comments,
                     conditional_comments),
            by = 'tag_code') %>%
  #filter out returning adults by limiting observations to the calendar year following tagging
  filter(first_obs <= paste0(tag_yr+1, '-12-31')) 

ggplot(ch_long, aes(x = first_obs)) +
  geom_histogram() +
  facet_wrap(~release_site, scale = "free_y")

#save year-specific ch for checking and summaries

saveRDS(ch_long, file = paste0("output/estimates/ch_parr_", tag_yr, ".RDS"))


ch_long_list <- ch_long %>%
  group_by(srr, release_site) %>%
  nest() %>%
  left_join(obs_tbl)



#define all as summer chinook so they will run together
ch_long_list <- ch_long_list %>%
  mutate(srr = '12W')

ch_long_list <- ch_long_list %>%
  mutate(
    ch = map2(
      data,
      obs_loc,
      ~ build_capture_histories(
        tag_history = .x,
        locs_def = .y,
        site_col = "obs_site",
        tag_col = "tag_code",
        time_col = "first_obs",
        censor_col = "last_censored",
        covariate_cols = NULL
      )
    )
  )

marray_list <- ch_long_list %>%
  mutate(m_array = map(
    ch, ~(marray = .x$m_array))) %>%
  select(srr, release_site, m_array) %>%
  rename(release_group = release_site)



#full shared detection site model
shared_p <- "All"
shared_phi <- "None"

fit_marray <- fit_marray_cjs_multiple(marray_list, 
                                        shared_p, shared_phi)

saveRDS(fit_marray, file = paste0("output/estimates/parr", tag_yr, ".RDS"))

#fit_ind <- run_cjs_model_set(ch_long_list %>% mutate(release_season = "summer"))

