# Purpose: Load mark and interogation files to develop a survival estimate.
# Author: Ryan N. Kinzer
# Created Date: April 28, 2026


# load libraries
library(PITmodelR)
library(tidyverse)

source('./scripts/run_cjs_model_set.R')

# get interrogation site configuration
site_config <- DART_site_config()

# get Data ---- 

activity <- 'rst' # used for directory path and hatchery, rst
mig_yr = 2025

# check for mark groups and if interrogation files exist (requires two specific file names "....Mark.rds" and "....INT.rds")
file_df <- check_mark_groups(activity = activity, migratory_year = mig_yr)

# we can only run survivals if we have interrogation or complete tag history information
mark_groups <- file_df %>% filter(cth_exists)

no_cth <- file_df %>% filter(!cth_exists)

# Need to map across all the groups contained in mark_groups

mark_df <- map_df(mark_groups$mark_file,
                  ~readRDS(file.path(file_path, .))
                  )

spp <- c('11W', '12W', '32W', '25W')
mark_site <- c('IMNTRP', 'JOHTRP', 'LOLTRP', 'SECTRP')

# list of tag_ids to use as marks
tag_ids <- mark_df %>%
  filter(species_run_rear_type %in% spp,
         release_site %in% mark_site,
         !grepl('Y', text_comments),
         !grepl('Y|M', conditional_comments)) %>%
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

# need to check MRR and interrogation site metadata for appropriate use
#mrr_sites <-..........

#obs_sites <- get_interrogation_sites(active_only = FALSE)
# obs_sites <- PITcleanr::buildConfig()

# filter to INT sites in tag_history$site_code

# my_sites = obs_sites %>%
#   filter(site_code %in% unique(ch_long$obs_site)) %>%
#   select(contains('site')) %>%
#   distinct()

# some MRR sites are avian bird colonies
# we now need to change the recovery observations to censored to reflect fish
# being removed

# ch_long$last_censored[ch_long$last_disposition == 'Recovery'] <- TRUE

# what sites should we include in the analysis

lgr <- c('GRJ', 'GRS')
down <- c("GOJ", "LMJ", "ICH", "MCJ", "JDJ", "B2J", "BCC",
          "PD5", "PD6", "PD7", "PD8", "PDO", "PDW", "TWX")

obs_tbl <- tibble(
  release_site = c("IMNTRP", "JOHTRP", "LOLTRP", "SECTRP"),
  obs_loc = list(
    list(
      IMNTRP = "IMNTRP",
      LGR = lgr,
      Down = down
    ),
    list(
      JOHTRP = "JOHTRP",
      ESS = "ESS",
      SFG = "SFG",
      LGR = lgr,
      Down = down
    ),
    list(
      LOLTRP = "LOLTRP",
      LGR = lgr,
      Down = down
    ),
    list(
      SECTRP = "SECTRP",
      ZEN = "ZEN",
      SFG = "SFG",
      LGR = lgr,
      Down = down
    )
  )
)

ch_long_list <- ch_long %>%
  group_by(srr, release_site, release_season) %>%
  nest() %>%
  left_join(obs_tbl)

model_results <- run_cjs_model_set(ch_long_list)

final_lgr_df <- model_results$estimates
model_lgr_tbl <- model_results$model_tbl

saveRDS(final_df, paste0(file_path,'/survival_estimates_LGR.rds'))


#---------- Old stuff--single dataset for testing

obs_locs_johtrp <- obs_tbl$obs_loc[[2]]

ch_list <- build_capture_histories(
  tag_history = ch_long,
  locs_def = obs_locs_johtrp,
  site_col = "obs_site",
  tag_col = "tag_code",
  time_col = "first_obs",
  censor_col = "last_censored",
  covariate_cols = NULL #c("grp")
)

#View(ch_list$ch_freq)
#View(ch_list$m_array)


# the first model from marked assumes all fish are available for detection so we
# need to change removals to observed and released
cjs_data <- ch_list$ch_data
cjs_data$ch <- gsub('2', '1',cjs_data$ch)

# this will result in an error because I have censored observations.
# fit the model with non-censored CJS
fit_cjs <- fit_marked_cjs(cjs_data, #ch_list$ch_data,
                      phi_formula = ~ time,
                      p_formula   = ~ time,
                      hessian = TRUE,       # for uncertainty
                      conf_level  = 0.95)

fit_cjs$phi
fit_cjs$p
fit_cjs$cum_phi


# fit model using the multi-state model
ms_data <- build_multistate_histories(ch_list$ch_data)

fit_mscjs <- fit_marked_mscjs(
  ms_data = ms_data,
  s_formula = ~ time,
  p_formula = ~ time,
  psi_formula = ~ -1 + stratum:tostratum
)

fit_mscjs$phi
fit_mscjs$p
fit_mscjs$psi
fit_mscjs$cum_phi


# model cjs using the m-array from Burnham
fit_marray <- fit_marray_cjs(
  m_array = ch_list$m_array
)

fit_marray$phi
fit_marray$p
fit_marray$lambda
fit_marray$cum_phi









# Interrogation Detail File
int_df <- read.csv(file.path(file_path, mark_groups$int_file[1])) %>%
  select(pit_tag = Tag.Code,
         obs_site = Site.Code.Value,
         obs_name = Site.Info.Name,
         obs_datetime = Obs.Time.Value,
         antenna_id = Antenna.ID,
         obs_lat = Site.Latitude.Value,
         obs_long = Site.Longitude.Value,
         obs_rkm = Site.RKM.Value,
         contains('Tag.'),
         release_lat = Release.Site.Latitude.Value,
         release_long = Release.Site.Longitude.Value) %>%
  mutate(
    obs_datetime = mdy_hms(obs_datetime),
    obs_date = as.Date(obs_datetime)) %>%
  arrange(pit_tag, obs_date) %>%
  mutate(row_id = row_number()) %>%
  group_by(pit_tag) %>%
  mutate(
    event = cumsum(obs_site != lag(obs_site, default = first(obs_site))) + 1L
  ) %>%
  ungroup()


names(int_df) <- gsub("\\.", "_", tolower(names(int_df)))

int_df <- int_df %>%
  inner_join(site_config, by = c("obs_site", "antenna_id")) %>%
  filter(between(obs_date, date_start, date_end))

obs_df <- int_df %>%
  group_by(pit_tag, obs_site, event) %>%
  summarise(
    first_obs = min(obs_datetime),
    last_obs  = max(obs_datetime),
    last_disposition = last(disposition[order(obs_datetime)]),
    last_censored = last(censored[order(obs_datetime)]),
    .groups = "drop"
  )

# need to create a data frame that binds_rows of marks with int.

ch_long <- mark_df %>%
  select(pit_tag, obs_site = release_site, obs_date = release_date) %>%
  mutate(last_disposition = 'Returned to River',
         last_censored = FALSE) %>%
  bind_rows(obs_df %>%
              select(pit_tag, obs_site, obs_date = last_obs, last_disposition, last_censored)) %>%
  arrange(pit_tag, obs_date)


table(ch_long$obs_site)


# R. Vosbigian's functions from "space4time"






# read_DART_file <- function(filepath, ...) {
#   dots <- list(...)
#   if (length(dots)  == 0) {
#     files = list(filepath)
#   } else {
#     # stopifnot(is.character(...))
#     files = list(filepath,...)
#   }
#
#
#   # find first row:
#   skip_df = read.csv(files[[1]])
#
#   skip = grep("[#]RelGrpStartDate",skip_df[,1])
#
#   capture_data <- readr::read_csv(files[[1]],skip = skip) %>%
#     as.data.frame() %>%
#     dplyr::filter(!is.na(TagID))
#
#
#
#   if (length(files) > 1) {
#     for (i in 2:length(files)) {
#       # find first row:
#       skip_df = read.csv(files[[1]])
#
#       skip = grep("[#]RelGrpStartDate",skip_df[,1])
#
#       tmp <- readr::read_csv(filepath,comment = "##") %>%
#         as.data.frame() %>%
#         dplyr::filter(!is.na(TagID))
#
#       capture_data <- rbind(capture_data,tmp)
#
#     }
#   }
#
#   capture_data$time <- tryCatch(lubridate::as_datetime(capture_data$RelVTime),
#                                 warning = function(e){
#                                   tryCatch(lubridate::mdy(capture_data$RelVTime),
#                                            warning = function(e2){
#                                              lubridate::as_date(
#                                                stringr::str_split_i(capture_data$RelVTime,
#                                                                     " ",
#                                                                     1),
#                                                format = "%m/%d/%Y")
#                                            })
#
#
#                                 })
#
#
#   capture_data$RecTime <- tryCatch(lubridate::as_datetime(capture_data$ObsDateLast),
#                                 warning = function(e){
#                                   tryCatch(lubridate::mdy(capture_data$ObsDateLast),
#                                            warning = function(e2){
#                                              lubridate::as_date(
#                                                stringr::str_split_i(capture_data$ObsDateLast,
#                                                                     " ",
#                                                                     1),
#                                                format = "%m/%d/%Y")
#                                            })
#
#
#                                 })
#
#
#   initialreleases <- capture_data %>%
#     dplyr::mutate(#time = lubridate::as_datetime(RelVTime),
#            time = lubridate::year(time),
#            # RelSite = ifelse(RelSite == "BBCTRP","BIGBEC",RelSite),
#            removed = FALSE) %>%
#     dplyr::select(id = TagID,site = RelSite,time = time, removed = removed)  %>%
#     dplyr::distinct(id,site,time,removed)
#
#   recaps <- capture_data %>%
#     dplyr::filter(!is.na(ObsSite)) %>%
#     # dplyr::mutate(ObsSite = case_when(ObsSite %in% c("B2J","BCC") ~ "BON",
#     #                            TRUE ~ ObsSite)) %>%
#     dplyr::filter(TagID %in% initialreleases$id) %>%
#     dplyr::mutate(removed = FALSE) %>%
#     dplyr::mutate(time = RecTime,#lubridate::as_datetime(ObsDateLast),
#            time = lubridate::year(time)) %>%
#     dplyr::select(id = TagID,site = ObsSite,time = time, removed = removed)
#
#   combined_capture_data <- rbind(initialreleases,recaps) %>%
#     dplyr::arrange(id, site, time)
#
#   return(combined_capture_data)
# }


# fix no visible binding note
ObsDateLast <- ObsSite <- RelSite <- RelVTime <- TagID <- NULL


# source("C:/Users/rvosbigian/OneDrive - University of Idaho/Post_masters/cjs_s4t/process_site_info_01_06_2025.R")
identify_barged_fish <- function(capture_data,parsed_df) {
  
  N_nositeinfo <- 0
  
  # investigate removed individuals
  for (i in 1:nrow(capture_data)) {
    if (is.na(capture_data$ObsSite[i])) next()
    
    tmp_obssite <- capture_data$ObsSite[i]
    tmp_obsmonitor <- capture_data$ObsMonitor[i]
    tmp_time <- capture_data$time[i]
    
    if (!(tmp_obssite %in% c("GOJ","GRJ","LMJ"))) next()
    
    tmp_siteinfo <- parsed_df %>%
      dplyr::filter(sitecode == tmp_obssite,
                    grepl(tmp_obsmonitor,arrayname),
                    date_start < tmp_time,
                    date_end > tmp_time)
    
    if (nrow(tmp_siteinfo) == 0) {
      capture_data$removed[i] <- FALSE
      N_nositeinfo <- N_nositeinfo + 1
      # stop("1")
    }
    
    if (length(unique(tmp_siteinfo$releasecode1)) == 1) {
      capture_data$removed[i] <- ifelse(tmp_siteinfo$releasecode1[1] == 'T',TRUE,FALSE)
    } else {
      
      if (nrow(tmp_siteinfo) > 1) {
        tmp_siteinfo <- tmp_siteinfo %>%
          dplyr::filter(exception == TRUE)
        
        if (nrow(tmp_siteinfo) == 1) {
          capture_data$removed[i] <- ifelse(tmp_siteinfo$releasecode1 == 'T',TRUE,FALSE)
        } else {
          tmp_releasecode <-unique(tmp_siteinfo$releasecode1)
          if (length(tmp_releasecode) == 1) {
            capture_data$removed[i] <- ifelse(tmp_releasecode == 'T',TRUE,FALSE)
          } else {
            stop("2")
          }
          
        }
        
        
      }
      
      # capture_data$removed[i] <- ifelse(tmp_siteinfo$releasecode1 == 'T',TRUE,FALSE)
      
      
    }
    
  }
  
  if (N_nositeinfo > 0) {
    warning(paste0("Number of observations without info. Assuming not transported. N = ",N_nositeinfo))
  }
  
  
  return(capture_data$removed)
  
  
} # end function




#' Read DART Basin TribPit observation files and combine with auxiliary age data
#'
#' @description
#' Read DART Basin TribPit observation files and combine with auxiliary age data
#'
#'
#'
#'
#' @param filepath the path to the DART Basin TribPit observation file or a
#'     vector of filepaths.
#' @param aux_age_df A data frame object in the same format as aux_age_df
#'     (see documentation for `s4t_ch`)
#' @param DART_config the path to a site configuration file `"site_config.txt"` generated by
#'     Columbia Basin Research. If blank, defaults to `www.cbr.washington.edu/downloads/paramest/sites_config.txt`
#'
#' @returns A list with two elements. The first is the formatted capture history data (`ch_df`)
#'     and the second is the age auxiliary data (`aux_age_df`).
#'
#'
#' @details
#' Currently, only Dams in the Snake River are checked for transported fish.
#'
#'
#' @export
#'
#'
read_DART_file <- function(filepath,aux_age_df,DART_config = "https://www.cbr.washington.edu/downloads/paramest/sites_config.txt") {
  
  if (is(filepath,"list")) {
    files = filepath
  } else {
    # stopifnot(is.character(...))
    files = as.list(filepath)
  }
  
  
  age_df <- aux_age_df
  
  # find first row:
  skip_df = utils::read.csv(files[[1]])
  
  skip = grep("[#]RelGrpStartDate",skip_df[,1])
  
  suppressMessages(capture_data <- readr::read_csv(files[[1]],skip = skip) %>%
                     as.data.frame() %>%
                     dplyr::filter(!is.na(TagID))
  )
  
  
  if (length(files) > 1) {
    for (i in 2:length(files)) {
      # find first row:
      skip_df = utils::read.csv(files[[i]])
      
      skip = grep("[#]RelGrpStartDate",skip_df[,1])
      
      # skip_df = readLines(files[[i]])
      #
      # skip = grep("[#]RelGrpStartDate",skip_df)
      
      suppressWarnings(suppressMessages(tmp <- readr::read_csv(files[[i]],comment = "##",skip = skip) %>%
                                          as.data.frame() %>%
                                          dplyr::filter(!is.na(TagID))))
      
      capture_data <- rbind(capture_data,tmp)
      
    }
  }
  
  
  
  
  capture_data$time <- tryCatch(lubridate::as_datetime(capture_data$RelVTime),
                                warning = function(e){
                                  tryCatch(lubridate::mdy(capture_data$RelVTime),
                                           warning = function(e2){
                                             lubridate::as_date(
                                               stringr::str_split_i(capture_data$RelVTime,
                                                                    " ",
                                                                    1),
                                               format = "%m/%d/%Y")
                                           })
                                  
                                  
                                })
  
  
  
  
  
  capture_data$RecTime <- tryCatch(lubridate::as_datetime(capture_data$ObsDateLast),
                                   warning = function(e){
                                     tryCatch(lubridate::mdy(capture_data$ObsDateLast),
                                              warning = function(e2){
                                                lubridate::as_date(
                                                  stringr::str_split_i(capture_data$ObsDateLast,
                                                                       " ",
                                                                       1),
                                                  format = "%m/%d/%Y")
                                              })
                                     
                                     
                                   })
  
  # if (is.null(DART_config)) {
  #   DART_config <- "https://www.cbr.washington.edu/downloads/paramest/sites_config.txt"
  # }
  
  parsed_df <- process_site_config(DART_config)
  
  capture_data$removed <- FALSE
  capture_data$removed <- (identify_barged_fish(capture_data,parsed_df))
  
  initialreleases <- capture_data %>%
    dplyr::mutate(#time = lubridate::as_datetime(RelVTime),
      time = lubridate::year(time),
      # RelSite = ifelse(RelSite == "BBCTRP","BIGBEC",RelSite),
      removed = FALSE) %>%
    dplyr::select(id = TagID,site = RelSite,time = time, removed = removed)  %>%
    dplyr::distinct(id,site,time,removed)
  
  recaps <- capture_data %>%
    dplyr::filter(!is.na(ObsSite)) %>%
    # dplyr::mutate(ObsSite = case_when(ObsSite %in% c("B2J","BCC") ~ "BON",
    #                            TRUE ~ ObsSite)) %>%
    dplyr::filter(TagID %in% initialreleases$id) %>%
    # dplyr::mutate(removed = FALSE) %>%
    dplyr::mutate(time = RecTime,#lubridate::as_datetime(ObsDateLast),
                  time = lubridate::year(time)) %>%
    dplyr::select(id = TagID,site = ObsSite,time = time, removed = removed)
  
  ### keep the first observation (typically kelts make up the later observations)
  recaps <- recaps %>%
    dplyr::group_by(id, site) %>%
    dplyr::arrange(id, site, time) %>%
    dplyr::mutate(N = 1:dplyr::n()) %>%
    dplyr::filter(N == 1) %>%
    dplyr::select(id, site, time, removed)
  
  
  combined_capture_data <- rbind(initialreleases,recaps) %>%
    dplyr::arrange(id, site, time)
  
  
  
  # adding rows to age_df if missing
  missing_obs_aux <- setdiff(capture_data$TagID,age_df$id)
  if (length(missing_obs_aux) > 0) {
    needed_cols <- c("obs_time","obs_site","id","ageclass","FL")
    check_cols <- setdiff(needed_cols,colnames(age_df))
    
    if (length(check_cols) > 0) {
      message(paste0("These columns are missing from age_df:",
                     paste0(check_cols,collapse = ", ")))
      stop("aux_age_df does not contain necessary columns")
    }
    
    needed_cap_cols <- c("RelSite","TagID","Lgth")
    check_cap_cols <- setdiff(needed_cap_cols,colnames(capture_data))
    if (length(check_cap_cols) > 0) {
      message(paste0("These columns are missing from DART file:",
                     paste0(check_cols,collapse = ", ")))
      stop(paste0("DART fils do not contain necessary columns"))
    }
    
    row_na_df <- age_df[1,]
    row_na_df[,1:ncol(row_na_df)] <- NA
    if (length(missing_obs_aux) > 1) {
      row_na_df <- row_na_df %>% tibble::add_row(id = rep(NA,length(missing_obs_aux)-1))
    }
    
    
    na_capture_data <- capture_data %>%
      dplyr::filter(TagID %in% missing_obs_aux) %>%
      dplyr::distinct(TagID,RelSite,time,.keep_all = TRUE)
    
    
    row_na_df$obs_site <- na_capture_data$RelSite
    
    row_na_df$id <- na_capture_data$TagID
    
    row_na_df$obs_time <- lubridate::year(na_capture_data$time)
    
    row_na_df$FL <- na_capture_data$Lgth
    
    message(paste0("Appending N = ",length(missing_obs_aux)," observations to age_df. \n    Using RelSite as obs_site, year of release as obs_time, and Lgth as FL"))
    
    message(paste0("Note: the following columns were populated with NAs: \n    ",
                   paste0(setdiff(colnames(row_na_df),c("obs_site",
                                                        "id",
                                                        "obs_time",
                                                        "FL")),
                          collapse = ", "
                   )
    )
    )
    
    age_df <- rbind(age_df,row_na_df)
  }
  
  # capture_data
  
  
  return(list(ch_df = combined_capture_data,
              aux_age_df = age_df))
}


#' Remove observations of kelts
#'
#' @description
#' Removes any observations following an observation at specified sites. Note that
#'     the observations at the specified sites are retained.
#'
#'
#' @export
#'
#' @param ch_df a capture history data frame (see documentation for s4t_ch)
#' @param kelt_obssite a vector of sites that identify kelts (i.e. adult fish ladders)
#'
#' @returns a capture history data frame without observations following
#'     observations at the specified sites.
#'
#' @export
#'
#'
remove_kelt_obs <- function(ch_df, kelt_obssite) {
  kelt_sites <- intersect(as.character(kelt_obssite),unique(ch_df$site))
  
  message(paste0("Dropping observations after observations at the following sites: ",paste0(kelt_sites,collapse = ", ")))
  
  
  kelt_obs <- ch_df %>%
    dplyr::filter(site %in% kelt_sites) %>%
    dplyr::group_by(id) %>%
    dplyr::arrange(time) %>%
    dplyr::summarize(kelt_time = dplyr::first(time))
  
  remove_kelts_ch_df <- ch_df %>%
    dplyr::left_join(kelt_obs, by = "id") %>%
    dplyr::group_by(id) %>%
    dplyr::filter(is.na(kelt_time) | time <= kelt_time) %>%
    dplyr::ungroup() %>%
    dplyr::select(id, site, time, removed)
  
  return(as.data.frame(remove_kelts_ch_df))
  
}



estimate_cohort_surv_vc <- function(obj) {
  res <- obj$res
  s4t_ch <- obj$fit$s4t_ch
  obj$fit$indices_theta$g
  
  n_init_relsite <- s4t_ch$ch_info$n_init_relsite
  n_groups <- length(unique(obj$fit$indices_theta$g)) # may cause issues
  
  
  
  init_relsite_list <- s4t_ch$ch_info$init_relsite_list
  set_min_a <- s4t_ch$ch_info$observed_relative_min_max$set_min_a
  set_max_a <- s4t_ch$ch_info$observed_relative_min_max$set_max_a
  
  
  max_s_rel <- s4t_ch$ch_info$max_s_rel
  max_t_recap <- s4t_ch$ch_info$max_t_recap
  n_sites <- s4t_ch$ch_info$n_sites
  recap_sites <- s4t_ch$ch_info$recap_sites
  last_sites <- s4t_ch$ch_info$last_sites
  recap_sites_not_last <- s4t_ch$ch_info$recap_sites_not_last
  not_last_sites <- c(1:n_sites)[!(1:n_sites %in% s4t_ch$ch_info$last_sites)]
  sites_config <- s4t_ch$s4t_config$sites_config
  holdover_config <- s4t_ch$s4t_config$holdover_config
  
  
  
  mod_mat_theta <- obj$fit$mod_mat_theta
  indices_theta <- obj$fit$indices_theta
  
  min_ageclass_mat <- s4t_ch$ch_info$min_ageclass_mat
  max_ageclass_mat <- s4t_ch$ch_info$max_ageclass_mat
  
  ja <- numDeriv::jacobian(function(x, ...) stats::plogis(return_cohort_surv(par = x, ...)),
                           x = res$par,
                           n_init_relsite = n_init_relsite,
                           init_relsite_list = init_relsite_list,
                           n_groups = n_groups,
                           set_min_a = set_min_a,
                           set_max_a = set_max_a,
                           max_s_rel = max_s_rel,
                           max_t_recap = max_t_recap,
                           n_sites = n_sites,
                           recap_sites = recap_sites,
                           last_sites = last_sites,
                           recap_sites_not_last = recap_sites_not_last,
                           not_last_sites = not_last_sites,
                           sites_config = sites_config,
                           holdover_config = holdover_config,
                           mod_mat_theta = mod_mat_theta,
                           indices_theta = indices_theta,
                           min_ageclass_mat = min_ageclass_mat,
                           max_ageclass_mat = max_ageclass_mat)
  
  vc <- solve(res$hessian)
  
  der_vc <- ja %*% vc %*% t(ja)
  
  return(der_vc)
  
}




#' Calculate abundance estimates
#'
#' @description
#' Calculate abundance estimates of cohorts
#'
#'
#' @export
#'
#' @param obj a capture history data frame (see documentation for s4t_ch)
#' @param abund a vector of sites that identify kelts (i.e. adult fish ladders)
#' @param type the type of summarization for the abundance estimates. BroodYear
#'     summarizes by broodyear (s - a1), ReleaseYear summarizes by time (s),
#'     and None does not summarize.
#'
#' @returns a capture history data frame without observations following
#'     observations at the specified sites.
#'
#' @export
#'
#'
abundance_estimates <- function(obj, abund, type = c("BroodYear","ReleaseYear","None")) {
  
  abund <- as.data.frame(abund)
  
  tmp_mis <- setdiff(c("a1","s","j","abundance","abundance_se"),colnames(abund))
  if ("abundance_se" %in% tmp_mis) {
    abund$abundance_se <- 0
  }
  
  tmp_mis <- setdiff(c("a1","s","j","abundance","abundance_se"),colnames(abund))
  
  if (length(tmp_mis) > 0) {
    stop(paste0("Missing columns from `abund`: ",paste0(tmp_mis,collapse = ", ")))
  }
  
  tmp_add <- setdiff(colnames(abund),c("a1","s","j","r","g","abundance","abundance_se"))
  
  if (length(tmp_add) > 0) {
    message(paste0("Only the following columns can be included in `abund`: ",
                   paste0(c("a1","s","j","r","g","abundance","abundance_se"),collapse = ", ")))
    stop(paste0("Remove columns from `abund`: ",paste0(tmp_add,collapse = ", ")))
  }
  
  abund$j <- as.character(abund$j)
  
  
  
  cohort_transitions <- obj$cohort_transitions
  suppressMessages(transition_df <- dplyr::left_join(cohort_transitions,abund))
  
  if (is(obj,"s4t_cjs_ml")) {
    
    
    
    # cohort_vc <- obj$fit$cohort_surv_vc
    cohort_vc <- estimate_cohort_surv_vc(obj = obj)
    
    
    
  } else if (is(obj,"s4t_cjs_rstan")) {
    colnames(transition_df)[colnames(transition_df) == "mean"] <- "estimate"
    colnames(transition_df)[colnames(transition_df) == "mean"] <- "estimate_se"
    
    tmp_cohort_trans_parnames <- paste0("cohort_surv[",1:obj$res@sim$dims_oi$cohort_surv,"]")
    # rstan::as.matrix.stanfit(obj$res,pars = tmp_cohort_trans_parnames)
    tmp_cohort <- rstan::extract(obj$res,pars = tmp_cohort_trans_parnames)
    cohort_vc <- stats::cov(do.call(cbind,tmp_cohort))
    
    # stats::cov()
  }
  
  
  
  # transition_df$abundance_se <- runif(nrow(transition_df),20,50)
  
  abund_vc <- diag(transition_df$abundance_se^2)
  
  
  
  outer_N <- transition_df$abundance %*% t(transition_df$abundance)
  outer_Theta <- transition_df$estimate %*% t(transition_df$estimate)
  
  cohort_abund_vc <- (abund_vc * cohort_vc) +
    (abund_vc * outer_Theta) +
    (outer_N * cohort_vc)
  
  # for approximation, just set NaNs to zero:
  if (sum(is.nan(cohort_abund_vc) > 0)) {
    message("Setting Nan values in vcov matrix to 0, results are approximate.")
    cohort_abund_vc <- ifelse(is.nan(cohort_abund_vc),0,cohort_abund_vc)
  }
  if (sum(is.na(cohort_abund_vc) > 0)) {
    message("Setting Nan values in vcov matrix to 0, results are approximate.")
    cohort_abund_vc <- ifelse(is.na(cohort_abund_vc),0,cohort_abund_vc)
  }
  
  
  calc_se <- function(Index,cohort_abund_vc) {
    
    a = matrix(rep(0,nrow(cohort_abund_vc)),nrow=1)
    a[Index] <- 1
    as.numeric(sqrt(a %*% cohort_abund_vc %*% t(a)))
  }
  
  if (type == "BroodYear") {
    
    
    
    
    
    transition_df %>%
      dplyr::mutate(Index = 1:dplyr::n(),
                    broodyear = s - a1) %>%
      dplyr::mutate(cohort_abund = estimate * abundance,
                    cohort_abund_se = sqrt((abundance_se^2 + abundance^2) * (diag(cohort_vc) + estimate^2) - ((abundance^2) * (estimate^2)))) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(broodyear,r,g,j) %>% #View()
      dplyr::summarize(abundance_broodyear = sum(cohort_abund),
                       # paste0(Index,collapse = " "),
                       abundance_broodyear_se = calc_se(Index,cohort_abund_vc)
      )
    
    
  } else if (type == "None") {
    transition_df %>%
      dplyr::mutate(Index = 1:dplyr::n(),
                    broodyear = s - a1) %>%
      dplyr::mutate(estimate_se = sqrt(diag(cohort_vc)),
                    cohort_abund = estimate * abundance,
                    cohort_abund_se = sqrt((abundance_se^2 + abundance^2) * (diag(cohort_vc) + estimate^2) - ((abundance^2) * (estimate^2))))
  } else if (type == "ReleaseYear") {
    transition_df %>%
      dplyr::mutate(Index = 1:dplyr::n(),
                    broodyear = s - a1) %>%
      dplyr::mutate(cohort_abund = estimate * abundance,
                    cohort_abund_se = sqrt((abundance_se^2 + abundance^2) * (diag(cohort_vc) + estimate^2) - ((abundance^2) * (estimate^2)))) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(s,j,r,g) %>%
      dplyr::summarize(abundance_releaseyear = sum(cohort_abund),
                       abundance_releaseyear_se = calc_se(Index,cohort_abund_vc))
    
  }
  
}