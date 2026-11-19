#Retrieve mark files for RST and parr tagging for comparison to 2025 parr tagging estimates

# load libraries
#remotes::install_github("ryankinzer/PITmodelR", ref = "dev", build_vignettes = F, force = T)
library(PITmodelR)
library(tidyverse)

# parameters

#RST data files

tag_yr <- c(2025, 2024, 2023, 2022, 2021) # could be multiple years c(2025, 2026)

proj_codes <- c(
  'CDR', 
  'OGR', 
  'SFC'
 ) 

proj_years <- expand.grid('mrr_project' = proj_codes,
                          'year' = tag_yr,
                          KEEP.OUT.ATTRS = FALSE,
                          stringsAsFactors = FALSE)

#think about adding specific srr codes here

# retrieve mrr file info for each project code and year

safe_get_mrr <- possibly(
  ~ get_mrr_files(code = .x, year = .y),
  otherwise = NULL
)

all_files_ls <- map2(
  proj_years$mrr_project,
  proj_years$year,
  function(x, y) {
    out <- safe_get_mrr(x, y)
    if (is.null(out)) warning("Warning: Failed for project code: ", x, " and year: ", y)
    out
  }
) %>%
  compact()

all_files <- all_files_ls %>%
  bind_rows() %>%
  mutate(user_ext = str_extract(name, "(?<=-)([A-Za-z0-9]{3})(?=\\.)"))

# check file types
table(all_files$projectCode, all_files$fileTypeExtension)
table(all_files$user_ext, all_files$projectCode)


# get mrr file data

mrr_dat_ls <- get_batch_file_data(all_files$name,
                                  pdvs = "attach",
                                  map_pdvs = FALSE)

saveRDS(mrr_dat_ls, file = paste0('./data/output/mrr_files/johtrp_mrr_files_',gsub(', ',"_",toString(tag_yr)),'.rds'))

### Parr tagging files

tag_yr <- c(2025, 2024, 2023, 2022, 2021) # could be multiple years c(2025, 2026)

proj_codes <- c("GAA", "OGR")

proj_years <- expand.grid('mrr_project' = proj_codes,
                          'year' = tag_yr,
                          KEEP.OUT.ATTRS = FALSE,
                          stringsAsFactors = FALSE)

#think about adding specific srr codes here

# retrieve mrr file info for each project code and year

safe_get_mrr <- possibly(
  ~ get_mrr_files(code = .x, year = .y),
  otherwise = NULL
)

all_files_ls <- map2(
  proj_years$mrr_project,
  proj_years$year,
  function(x, y) {
    out <- safe_get_mrr(x, y)
    if (is.null(out)) warning("Warning: Failed for project code: ", x, " and year: ", y)
    out
  }
) %>%
  compact()

all_files <- all_files_ls %>%
  bind_rows() %>%
  mutate(user_ext = str_extract(name, "(?<=-)([A-Za-z0-9]{3})(?=\\.)"))

# check file types
table(all_files$projectCode, all_files$fileTypeExtension)
table(all_files$user_ext, all_files$projectCode)


# get mrr file data

mrr_dat_ls <- get_batch_file_data(all_files$name,
                                  pdvs = "attach",
                                  map_pdvs = FALSE)

saveRDS(mrr_dat_ls, file = paste0('./data/output/mrr_files/GAA_OGR_mrr_files_',gsub(', ',"_",toString(tag_yr)),'.rds'))





# parse events for RST

tag_yr <- c(2025, 2024, 2023)

mrr_dat_ls <- readRDS(paste0('./data/output/mrr_files/rst_mrr_files_',gsub(', ',"_",toString(tag_yr)),'.rds'))

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

# subset parr tagging events

check_groups <- c('event_year', 'event_season', 'migration_year',
                  'mrr_project', 'file_ext', 'session_message',
                  'capture_method', 'event_site', 'release_site',
                  'species_run_rear_type', 'event_type')

#write function to separate events by RST, hatchery, and parr tagging? 
events_rst <- mrr_events %>%
  #filter(migration_year == mig_yr) %>%
  filter(capture_method %in% c("SCREWT"),
         event_type != 'Passive Recapture')

rst_qa <- check_events(events_rst, check_groups)

rst_pdv_map <- inner_join(pdv_map, events_rst %>%
                             select(file_name) %>%
                             distinct()) %>%
  mutate(files =paste0(str_sub(file_name, 1, 3),'-',str_extract(file_name, "(?<=-)([A-Za-z0-9]{3})(?=\\.)")))


label_files_diff_rst <- rst_pdv_map %>%
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

definition_files_diff_rst <- rst_pdv_map %>%
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

#check date range
max(events_rst$event_date)
min(events_rst$event_date)

#check release sites
table(events_rst$release_site)

 
# Subset events to 'Mark' only, save output, and create tag id txt file for
# PTAGIS upload.  You could iterate across the event_type_value to create
# separate files for each event.

tag_yr <- 2021

#grouped by release site
#filter for wild sp/sum Chinook
 save_event_groups(
   events = events_rst %>% filter(species_run_rear_type %in% c("12W", "11W")),
   event_type_value = "Mark",
   group_cols = c("event_year", "mrr_project", "release_site"),
   save_TagIds = TRUE,
   out_dir = paste0("./output/events_rst/tag_yr",tag_yr)
 )

 #Pull CTH files from PTAGIS using tagid text files
 #add Antenna, Event Date Time, and Event Release Date Time to PTAGIS query
 
 
 ############# Parse events for parr tagging ###############
 
 
 tag_yr <- c(2025, 2024, 2023, 2022, 2021)
 
 mrr_dat_ls <- readRDS(paste0('./data/output/mrr_files/GAA_mrr_files_',gsub(', ',"_",toString(tag_yr)),'.rds'))
 
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
 
 # subset parr tagging events
 
 check_groups <- c('event_year', 'event_season', 'migration_year',
                   'mrr_project', 'file_ext', 'session_message',
                   'capture_method', 'event_site', 'release_site',
                   'species_run_rear_type', 'event_type')
 
 #write function to separate events by RST, hatchery, and parr tagging? 
 events_parr <- mrr_events %>%
   #filter(migration_year == mig_yr) %>%
   filter(capture_method %in% c("BSEINE", "SHOCK"),
          event_type != 'Passive Recapture')
 
 parr_qa <- check_events(events_parr, check_groups)
 
 parr_pdv_map <- inner_join(pdv_map, events_parr %>%
                             select(file_name) %>%
                             distinct()) %>%
   mutate(files =paste0(str_sub(file_name, 1, 3),'-',str_extract(file_name, "(?<=-)([A-Za-z0-9]{3})(?=\\.)")))
 
 
 label_files_diff_parr <- parr_pdv_map %>%
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
 
 definition_files_diff_rst <- parr_pdv_map %>%
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
 
 #check date range
 max(events_parr$event_date)
 min(events_parr$event_date)
 
 #check release sites
 table(events_parr$release_site)
 
 
 # Subset events to 'Mark' only, save output, and create tag id txt file for
 # PTAGIS upload.  You could iterate across the event_type_value to create
 # separate files for each event.
 
 #2025 files are already compiled
 
 
 years <- seq(2021, 2024, by = 1)
 
 for (y in 1:length(years)) {
   
  tag_yr <- years[y]
 
 #grouped by release site
 #filter for wild sp/sum Chinook
 #separate by tag year
 save_event_groups(
   events = events_parr %>% filter(species_run_rear_type %in% c("12W", "11W"), 
                                   event_year == tag_yr),
   event_type_value = "Mark",
   group_cols = c("event_year", "mrr_project", "release_site"),
   save_TagIds = TRUE,
   out_dir = paste0("./output/events_parr/tag_yr",tag_yr)
 )
 }
 
 #Pull CTH files from PTAGIS using tagid text files
 #add Antenna, Event Date Time, and Event Release Date Time to PTAGIS query
 
 
#add 2025 OGR tag files
 
 save_event_groups(
   events = events_parr %>% filter(species_run_rear_type %in% c("12W", "11W"), 
                                   event_year == 2025, 
                                   mrr_project == "OGR"),
   event_type_value = "Mark",
   group_cols = c("event_year", "mrr_project", "release_site"),
   save_TagIds = TRUE,
   out_dir = paste0("./output/events_parr/tag_yr",2025)
 )
 
 