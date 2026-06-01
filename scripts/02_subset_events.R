# Purpose: Subset into hatchery, RST, and parr tagging eveents and output TagIds
# Author: Ryan N. Kinzer
# Created Date: April 16, 2025


# load libraries
library(PITmodelR)
library(tidyverse)

# Get Data ---- 

tag_yr = c(2024,2025)
mig_yr = 2025

mrr_dat_ls <- readRDS(paste0('./output/mrr_files/npt_mrr_files_',gsub(', ',"_",toString(tag_yr)),'.rds'))

# check for count/tally fields in PDV fields
pdv_map <- mrr_dat_ls$pdv_map
tally_fields <- pdv_map[which(grepl('count|tally|tallied|nfish|n_fish|number', tolower(pdv_map$definition))),]

tally_fields %>% 
  mutate(mrr_project = substr(file_name,1,3)) %>% 
    distinct(mrr_project, pdv_column)

mrr_events <- mrr_dat_ls$sessions %>%
  left_join(mrr_dat_ls$events,
            by = 'file_name') %>%
  mutate(tag_code = ifelse(pit_tag == '..........',NA,pit_tag),
         file_ext = str_extract(file_name, "(?<=-)([A-Za-z0-9]{3})(?=\\.)"),
         event_ym = floor_date(event_date, "month"),
         release_ym = floor_date(release_date, "month"),
         event_year = year(event_date),
         event_season = case_when(
           between(month(event_date), 1,6) ~ 'Spring',
           between(month(event_date), 7, 12) ~ 'Summer/Fall')
         )

# helper function to QA/QC

# if (nrow(tmp) > 0) {
#   stop("Wild origin fish were found in the hatchery dataset.")
# } else {
#   message("QA/QC Check - Passed")
# }

check_events <- function(.data, group_vars){
  
  safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
  safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  
  tmp <- .data %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n_files    = n_distinct(file_name),
      n_events   = n(),
      length_min = safe_min(length),
      length_max = safe_max(length),
      weight_min = safe_min(weight),
      weight_max = safe_max(weight),
      .groups = 'drop'
    ) 
}

# subset hatchery marking events: DIPNET - should always be hatchery marking events
events_hat <- mrr_events %>%
  filter(migration_year == mig_yr) %>%
  filter(capture_method == 'DIPNET') %>%
  filter(mrr_project != 'GAA')

check_groups <- c('event_year', 'event_season', 'migration_year',
             'mrr_project', 'file_ext', 'session_message',
             'capture_method', 'event_site', 'release_site',
             'species_run_rear_type', 'event_type')

hat_qa <- check_events(events_hat, check_groups)

hat_pdv_map <- inner_join(pdv_map, events_hat %>%
                            select(file_name) %>%
                            distinct()) %>%
  mutate(files =paste0(str_sub(file_name, 1, 3),'-',str_extract(file_name, "(?<=-)([A-Za-z0-9]{3})(?=\\.)")))


label_files_diff <- hat_pdv_map %>%
  group_by(level, pdv_column) %>%
  filter(n_distinct(label) > 1) %>%
  ungroup() %>%
  group_by(level, pdv_column, label) %>%
  summarise(
    files = paste(sort(unique(files)), collapse = "; "),
    n_files = n_distinct(files),
    .groups = "drop"
  ) %>%
  arrange(level, pdv_column, label)

definition_files_diff <- hat_pdv_map %>%
  group_by(level, pdv_column) %>%
  filter(n_distinct(definition) > 1) %>%
  ungroup() %>%
  group_by(level, pdv_column, definition) %>%
  summarise(
    files = paste(sort(unique(files)), collapse = "; "),
    n_files = n_distinct(files),
    .groups = "drop"
  ) %>%
  arrange(level, pdv_column, definition)


# subset natural marking events into RST and parr collections
events_rst <- mrr_events %>%
  filter(migration_year == mig_yr) %>%
  filter(capture_method %in% c("SCREWT"),
         event_type != 'Passive Recapture') %>%
  mutate(n_fish = as.numeric(pdv1))

rst_qa <- check_events(events_rst, check_groups)

rst_pdv_map <- inner_join(pdv_map, events_rst %>%
                            select(file_name) %>%
                            distinct()) %>%
  mutate(files =paste0(str_sub(file_name, 1, 3),'-',str_extract(file_name, "(?<=-)([A-Za-z0-9]{3})(?=\\.)")))


label_files_diff <- rst_pdv_map %>%
  group_by(level, pdv_column) %>%
  filter(n_distinct(label) > 1) %>%
  ungroup() %>%
  group_by(level, pdv_column, label) %>%
  summarise(
    files = paste(sort(unique(files)), collapse = "; "),
    n_files = n_distinct(files),
    .groups = "drop"
  ) %>%
  arrange(level, pdv_column, label)

definition_files_diff <- rst_pdv_map %>%
  group_by(level, pdv_column) %>%
  filter(n_distinct(definition) > 1) %>%
  ungroup() %>%
  group_by(level, pdv_column, definition) %>%
  summarise(
    files = paste(sort(unique(files)), collapse = "; "),
    n_files = n_distinct(files),
    .groups = "drop"
  ) %>%
  arrange(level, pdv_column, definition)

events_parr <- mrr_events %>%
  #filter(migration_year == mig_yr) %>%
  filter(capture_method %in% c("BSEINE", "SHOCK", "HOOK", "PSEINE"),
         event_type != 'Passive Recapture') %>%
  mutate(n_fish = as.numeric(pdv1))

parr_qa <- check_events(events_parr, check_groups)

parr_pdv_map <- inner_join(pdv_map, events_rst %>%
                            select(file_name) %>%
                            distinct()) %>%
  mutate(files =paste0(str_sub(file_name, 1, 3),'-',str_extract(file_name, "(?<=-)([A-Za-z0-9]{3})(?=\\.)")))


label_files_diff <- parr_pdv_map %>%
  group_by(level, pdv_column) %>%
  filter(n_distinct(label) > 1) %>%
  ungroup() %>%
  group_by(level, pdv_column, label) %>%
  summarise(
    files = paste(sort(unique(files)), collapse = "; "),
    n_files = n_distinct(files),
    .groups = "drop"
  ) %>%
  arrange(level, pdv_column, label)

definition_files_diff <- parr_pdv_map %>%
  group_by(level, pdv_column) %>%
  filter(n_distinct(definition) > 1) %>%
  ungroup() %>%
  group_by(level, pdv_column, definition) %>%
  summarise(
    files = paste(sort(unique(files)), collapse = "; "),
    n_files = n_distinct(files),
    .groups = "drop"
  ) %>%
  arrange(level, pdv_column, definition)

# Subset events to 'Mark' only, save output, and create tag id txt file for
# PTAGIS upload.  You could iterate across the event_type_value to create
# separate files for each event.

save_event_groups(
  events = events_hat,
  event_type_value = "Mark",
  group_cols = c("migration_year", "species_run_rear_type", "mrr_project", "release_site"),
  save_TagIds = TRUE,
  out_dir = paste0("./output/events_hatchery/MY",mig_yr)
)

save_event_groups(
  events = events_rst,
  event_type_value = "Mark",
  group_cols = c("migration_year", "mrr_project"),
  save_TagIds = TRUE,
  out_dir = paste0("./output/events_rst/MY",mig_yr)
)

map(.x = c('2025', '2026'),
    ~save_event_groups(
      events = events_parr %>% filter(migration_year == .x),
      event_type_value = "Mark",
      group_cols = c("migration_year"),
      save_TagIds = TRUE,
      out_dir = paste0("./output/events_parr/MY",.x)
    )
    )





