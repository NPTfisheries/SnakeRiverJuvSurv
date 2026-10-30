# Purpose: develop models to estimate survival for multiple release groups
# Author: Michelle A. Briggs
# Created Date: 18 August 2026

remotes::install_github("ryankinzer/PITmodelR", ref = "dev", build_vignettes = F, force = T)

# load libraries
library(PITmodelR)
library(tidyverse)
library(marked)
library(RTMB)

plot_theme <- readRDS("templates/plot_theme_v2.RDS")

source('./scripts/run_cjs_model_set.R')

# get interrogation site configuration
site_config <- DART_site_config()

# get Data ---- 

activity <- 'parr' # used for directory path and hatchery, rst
tag_yr <-  2025

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

#spp <- c('11W', '12W', '32W', '25W')
#mark_site <- c('IMNTRP', 'JOHTRP', 'LOLTRP', 'SECTRP')


#testing multiple release group models with summer 2025 Chinook parr tagging

spp <- c('11W', '12W')
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
         release_site %in% mark_site,
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
table(ch_long$text_comments) #check on what all the text comments are!!
table(ch_long$conditional_comments) 


lgr <- c('GRJ', 'GRS')
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

ch_long %>% left_join(obs_tbl, by = "release_site")


  


ch_long %>%
  group_by(srr, release_site, release_season) %>% 
  count()  

#in this case, code all sites with 12W srr so they all run together in a single model
#otherwise, 11w and 12w run separately
ch_long <- ch_long %>%
  mutate(srr = '12W') %>%
  filter(release_site != "BEAR6C")

ch_long_list <- ch_long %>%
  group_by(srr, release_season, release_site) %>%
  nest() %>%
  left_join(obs_tbl, by = "release_site")


#select smaller group of release groups that share SFG and LGR detection sites

ch_long_list_sf <- ch_long_list %>%
  filter(release_site %in% c(
                             "SUMTIC", 
                             "GROUSC", 
                             "JOHNSC", 
                             "BURNLC", 
                             "SAEFSF", 
                             "SECESR", 
                             "LAKEC", 
                             "SALRSF"
                             ))


results_mult <- run_cjs_model_set_multiple(ch_long_list_sf)
#running MS CJS with multiple release groups is slow with 8 release groups (5 - 10 min)

group_est <- results_mult$estimates %>%
  select(release_group, interval, metric, from_site, to_site, 
         cjs_est, cjs_lcl, cjs_ucl, 
         mscjs_est, mscjs_lcl, mscjs_ucl, 
         marray_est, marray_lcl, marray_ucl) %>%
  mutate(method = "grouped")


#check on RTMB function to make sure produces warning for singular fit, indicating model didn't converge
#usually this doesn't produce SE estimates, so probably would be evident anyways

model_results_ind <- run_cjs_model_set(ch_long_list_sf)



ind_est <- model_results_ind$estimates %>%
  rename(release_group = release_site) %>% 
  select(release_group, interval, metric, from_site, to_site, 
         cjs_est, cjs_lcl, cjs_ucl,
         mscjs_est, mscjs_lcl, mscjs_ucl, 
         marray_est, marray_lcl, marray_ucl) %>%
  mutate(method = "independent")
  
 method_compare <- rbind(group_est, ind_est) %>%
   pivot_longer(
     cols = matches("^(cjs|mscjs|marray)_(est|lcl|ucl)$"),
     names_to = c("model", ".value"),
     names_pattern = "(cjs|mscjs|marray)_(est|lcl|ucl)")
  

#compare estimates from independent and grouped model

survival_plot <- method_compare %>%
  filter(metric == "survival", to_site != "Down") %>%
  unite(col = "int_label", from_site, to_site, sep = " to ", remove = FALSE) %>%
  mutate(interval = as.numeric(interval)) %>%
  mutate(int_label = fct_reorder(int_label, interval)) %>%
  ggplot(aes(x = int_label, y = est, group = interaction(model, method))) + 
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.1, position = position_dodge(width = 0.3)) +
  geom_point(size = 4, 
             aes(fill = method, shape = model), 
             color = "black", 
             position = position_dodge(width = 0.3)) + 
  scale_shape_manual(values = c(21, 22, 23)) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) + 
  plot_theme + 
  facet_wrap(~release_group, scales = "free_x", ncol = 3) +
  scale_y_continuous(limits = c(0, 1)) + 
  theme(strip.background = element_rect(fill = "grey", color = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab("Survival") + 
  xlab(" ")
survival_plot

#still thinking about the best way to visualize this
detection_plot <- method_compare %>%
  filter(metric == "detection", from_site != "Down", !is.na(est)) %>%
  mutate(type = if_else(method == "grouped", "grouped", release_group)) %>%
  ggplot(aes(x = type, y = est, group = interaction(model, type))) +
  #geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.1, 
  #              position = position_jitterdodge(dodge.width = 0.2, jitter.width = 0.2)) +
  geom_pointrange(size = 1, aes(fill = method, shape = model, ymin = lcl, ymax = ucl), color = "black", 
             #position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.6)
             position = position_dodge(width = 0.6)) + 
  scale_shape_manual(values = c(21, 22, 23)) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) + 
  plot_theme + 
  facet_wrap(~from_site) +
  scale_y_continuous(limits = c(0, 1)) + 
  theme(strip.background = element_rect(fill = "grey", color = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab("Detection") + 
  xlab("Release group model")

detection_plot

#detection estimates become more precise
#survival estimates typically become more precise but not every time
#tested with set of parr tagging groups from SF Salmon, all had shared detections at SFG



#look at timing of passage by detection sites

#probably a really inefficient way to do this
site_passage <- data.frame()

for (i in 1:nrow(ch_long_list)) {
  release_group <- ch_long_list$release_site[[i]]
  dat <- ch_long_list$data[[i]]
  
  dat <- dat %>% select(obs_site, first_obs) %>%
    mutate(release_group = release_group) 
  #%>%
  #  mutate(obs_site = if_else(obs_site %in% lgr, "LGR", 
  #                            if_else(obs_site %in% down, "Down", obs_site)))
  
  site_passage <- rbind(site_passage, dat)
  
}

site_passage %>%
  filter(obs_site %in% unique(unlist(obs_tbl$obs_loc))) %>%
  mutate(obs_site = if_else(obs_site %in% lgr, "LGR", 
                            if_else(obs_site %in% down, "Down", obs_site))) %>%
  ggplot(aes(x = release_group, y = first_obs)) + 
  geom_boxplot() + 
  facet_wrap(~obs_site)








##### Testing
# with a single set of release groups

#this creates capture histories by srr, release site, and season
ch_long_list_sf <- ch_long_list_sf %>%
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


#single df with release site, occasion, and location
site_mapping <- ch_long_list_sf %>%
  mutate(site_map = map(ch, "mapping")) %>%
  select(release_site, site_map) %>%
  rename(release_group = release_site) %>%
  unnest(site_map) %>%
  ungroup() %>%
  select(release_group, occasion, occ_idx) %>%
  distinct() 

max_occ <- max(site_mapping$occ_idx)

#adjust occasion indexing for release groups with a lower number of occasions
# this is important for sharing information for sites downstream (i.e., LGR and Down)
site_mapping <- site_mapping %>%
  group_by(release_group) %>%
  mutate(occ_idx = occ_idx + (max_occ - max(occ_idx))) %>%
  ungroup() %>%
  rename(occ = occ_idx, 
         detect_site = occasion)

marked_ch_list <- ch_long_list_sf %>%
  group_by(srr, release_season) %>% 
  #combine release sites together, keep separate by srr and release season
  summarise(ch_data = list(map2_dfr(ch, release_site, ~ mutate(.x$ch_data, release_site = .y))), 
            .groups = "drop") %>%
  #determine length of longest ch by group
  mutate(ch_data = map(ch_data, ~ {max_len <- max(nchar(.x$ch)) 
  .x %>%
    #add leading 0s to shorter ch, all release groups end with detections at LGR and downstream
    mutate(ch = str_pad(ch, width = max_len, side = "left", pad = "0"))})
  )

#CJS
#select a single srr/season
#marked_ch <- marked_ch_list$ch_data[[4]]
marked_ch <- marked_ch_list$ch_data[[1]]
marked_ch$ch <- gsub('2', '1',marked_ch$ch)
marked_ch <- marked_ch %>%
  mutate(release_group = as.factor(release_site))
marked_ch <- data.frame(marked_ch)

fit_mult <- fit_marked_cjs_multiple(marked_ch, site_mapping = site_mapping)


#MS CJS

marked_ch <- build_multistate_histories(marked_ch)

proc <- marked::process.data(
  marked_ch,
  model = "hmmMSCJS",
  groups = "release_group",
  strata.labels = c("A", "C")
)

ddl <- marked::make.design.data(proc)

ddl$p <- dplyr::left_join(ddl$p, site_mapping, by = c("release_group", "occ"))
ddl$p <- ddl$p %>%
  dplyr::mutate(detect_site = if_else(is.na(detect_site), release_group, detect_site))

ddl$p$detect_site <- as.factor(ddl$p$detect_site)

ddl$p$fix <- NA
ddl$p$fix[ddl$p$stratum == "C" & ddl$p$occ == max(ddl$p$occ)] <- 1

#fix p to 0 for stratum == C for detections above Down, because censoring doesn't occur
#model runs slightly faster
#check that this assumption holds

ddl$p$fix[ddl$p$stratum == "C" & ddl$p$occ < max(ddl$p$occ)] <- 0

ddl$Psi$fix <- NA
ddl$Psi$fix[ddl$Psi$stratum == "C" & ddl$Psi$tostratum == "A"] <- 0
ddl$Psi$fix[ddl$Psi$stratum == "C" & ddl$Psi$tostratum == "C"] <- 1

s_formula   = ~ time*release_group
p_formula   = ~ stratum*detect_site*time
psi_formula = ~ -1 + stratum:tostratum

mod <- marked::crm(
  proc,
  ddl,
  model = "hmmMSCJS",
  model.parameters = list(
    S   = list(formula = s_formula),
    p   = list(formula = p_formula),
    Psi = list(formula = psi_formula)
  ),
  hessian = TRUE)


ms_fit_mult <- fit_marked_mscjs_multiple(marked_ch, site_mapping)
