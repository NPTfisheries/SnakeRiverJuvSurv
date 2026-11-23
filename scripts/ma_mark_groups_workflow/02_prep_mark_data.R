# PREPARE MARK DATA
#
# Author(s): Mike Ackerman, Michelle Briggs, and Ryan Kinzer
#
# Date Created: November 19, 2026
#   Last Modified: November 23, 2026
#
# Purpose: Compile cached PTAGIS MRR data, prepare juvenile release events,
#   summarize available data, and save documented outputs for mark-group creation.
#
# NOTES: Includes Mark and Recapture events. Mark group selection and biological
#   exclusions are handled separately in 03_create_mark_groups.R.

# clear environment
rm(list = ls())

# load necessary packages
library(tidyverse)
library(here)
library(PITmodelR)

#------------------------------------
# SELECT PROJECT CODES AND FILE YEARS
project_codes = c("BDA", "CDR", "GAA", "MAP", "PJC", "OGR", "SCS", "SFC")
file_years    = 2021L:2025L
pdv_option    = "drop"  # "drop" or "map"

# directory where cached PTAGIS MRR data are located
cache_dir = here("data", "cache", "ptagis", "mrr")

# label run and directory for output
prep_id    = paste0("prep_", format(Sys.time(), "%Y%m%d_%H%M"))
output_dir = here("data", "processed", "mark_prep", prep_id)

#------------------------
# COMPILE CACHED MRR DATA

# identify cached MRR files
cache_files = crossing(
  project_code = project_codes,
  file_year    = file_years
) |>
  mutate(
    cache_path = file.path(cache_dir, project_code, paste0("mrr_", project_code, "_", file_year, ".rds"))
  )

# require all requested caches
missing_files = cache_files |>
  filter(!file.exists(cache_path))

if (nrow(missing_files) > 0) {
  stop("Missing MRR caches:\n", paste(missing_files$cache_path, collapse = "\n"), "\nRun 01_retrieve_mrr.R or revise the selected projects/years.")
}

# read cached datasets
mrr_ls = map(cache_files$cache_path, readRDS)

# verify cache structures
if (!all(map_lgl(
  mrr_ls,
  ~ all(c("sessions", "events", "pdv_map") %in% names(.x))
))) {
  stop("One or more MRR caches have an unexpected structure.")
}

# record input files and their contents
input_manifest = cache_files |>
  transmute(
    project_code,
    file_year,
    cache_file = file.path("data", "cache", "ptagis", "mrr", project_code, basename(cache_path)),
    md5        = unname(tools::md5sum(cache_path)),
    n_sessions = map_int(mrr_ls, ~ nrow(.x$sessions)),
    n_events   = map_int(mrr_ls, ~ nrow(.x$events))
  )

# compile sessions and events; retain source file year on each event
mrr = list(
  sessions = map_dfr(mrr_ls, "sessions"),
  events = map2_dfr(
    mrr_ls, cache_files$file_year,
    ~ mutate(.x$events, cache_year = .y)
  )
)

# map PDVs only when requested
if (pdv_option == "map") {
  
  mrr$pdv_map = map_dfr(mrr_ls, "pdv_map")
  mrr = PITmodelR::map_pdvs_to_cols(mrr)
  
} else if (pdv_option == "drop") {
  
  mrr$sessions = mrr$sessions |> select(-matches("^spdv[0-9]+$"))
  mrr$events = mrr$events |> select(-matches("^pdv[0-9]+$"))
  
}

rm(mrr_ls)

# verify file_name can identify sessions and match event records
if (anyDuplicated(mrr$sessions$file_name)) {
  stop("Duplicate file_name values found in MRR sessions.")
}

if (!all(mrr$events$file_name %in% mrr$sessions$file_name)) {
  stop("Some MRR events have no matching session.")
}

#------------------
# PREPARE MARK DATA

# event variables needed for mark-group creation and QC
event_vars = c(
  "file_name", "sequence_number", "cache_year", 
  "species_run_rear_type", "pit_tag", "capture_method", "mark_method", "life_stage", 
  "event_type", "event_site", "event_date", 
  "release_site", "release_date",
  "organization", "location_rkmext", 
  "length", "weight",
  "tagger", "conditional_comments"
)

# optional event variables, including any desired mapped PDVs
optional_vars = c(
  "text_comments", "migration_year", "mark_temperature",
  "release_temperature", "holding_temperature"
)

# retain juvenile Mark and Recapture events with valid PIT tags
candidate_marks = mrr$events |>
  filter(
    event_type %in% c("Mark", "Recapture"),
    life_stage == "Juvenile",
    !is.na(pit_tag),
    str_trim(pit_tag) != "",
    pit_tag != ".........."
  ) |>
  select(all_of(event_vars), any_of(optional_vars)) |>
  left_join(
    mrr$sessions |> select(file_name, mrr_project),
    by = "file_name",
    relationship = "many-to-one"
  ) |>
  mutate(
    species_run = case_when(
      str_sub(species_run_rear_type, 1, 2) %in% c("11", "12") ~ "SpSu_Chinook",
      str_sub(species_run_rear_type, 1, 2) == "13"            ~ "Fall_Chinook",
      str_sub(species_run_rear_type, 1, 1) == "3"             ~ "Steelhead",
      TRUE ~ "Other"
    ),
    rear = str_sub(species_run_rear_type, 3, 3),
    event_year = lubridate::year(event_date),
    release_year = lubridate::year(release_date),
    release_month = lubridate::month(release_date)
  ) |>
  filter(
    species_run != "Other",
    rear %in% c("H", "W")
  )

# original capture-method summary, before further classifications
candidate_mark_summary = candidate_marks |>
  count(
    species_run, rear, release_year, capture_method,
    name = "n_events"
  ) |>
  arrange(species_run, rear, release_year, desc(n_events))

#---------------------------------------
# DEFINE GROUPING VARIABLES AND QC FLAGS

# PTAGIS conditional-comment codes; see https://www.ptagis.org/Resources/ValidationCodes
mortality_codes = c("M", "MB", "MK", "MS", "SM")
yearling_codes  = c("Y", "PR")
concern_codes   = c("B", "I", "DB", "L")

# patterns match complete codes
mortality_pattern = paste0("\\b(", paste(mortality_codes, collapse = "|"), ")\\b")
yearling_pattern  = paste0("\\b(", paste(yearling_codes, collapse = "|"), ")\\b")
concern_pattern   = paste0("\\b(", paste(concern_codes, collapse = "|"), ")\\b")

candidate_marks = candidate_marks |>
  mutate(
    # capture-method groups
    method_group = case_when(
      capture_method %in% c("BSEINE", "PSEINE", "HOOK", "SHOCK") ~ "ROVING",
      capture_method == "SCREWT"                                 ~ "SCREWT",
      capture_method == "DIPNET"                                 ~ "DIPNET",
      TRUE ~ "OTHER"
    ),
    # species-specific seasons using release month
    season = case_when(
      species_run == "SpSu_Chinook" & between(release_month, 1, 6)  ~ "spring",
      species_run == "SpSu_Chinook" & between(release_month, 7, 8)  ~ "summer",
      species_run == "SpSu_Chinook" & between(release_month, 9, 12) ~ "fall",
      species_run == "Fall_Chinook" & between(release_month, 1, 12) ~ "subyearling",
      species_run == "Steelhead"    & between(release_month, 1, 6)  ~ "spring",
      species_run == "Steelhead"    & between(release_month, 7, 8)  ~ "summer",
      species_run == "Steelhead"    & between(release_month, 9, 12) ~ "fall",
      is.na(release_month) ~ "missing_date",
      TRUE ~ "outside_season"
    ),
    # biological concern flags; no exclusions applied yet
    flag_mortality = str_detect(str_to_upper(coalesce(conditional_comments, "")),  mortality_pattern),
    flag_yearling  = species_run == "SpSu_Chinook" & str_detect(str_to_upper(coalesce(conditional_comments, "")), yearling_pattern),
    flag_concern   = str_detect(str_to_upper(coalesce(conditional_comments, "")), concern_pattern)
  )

# remove compiled MRR tables
rm(mrr)

#------------------------------
# SUMMARIZE AVAILABLE MARK DATA

# availability before biological exclusions or minimum-size filtering
mark_availability = candidate_marks |>
  group_by(
    species_run, rear, release_year, method_group, release_site, season
  ) |>
  summarise(
    n_events    = n(),
    n_marks     = sum(event_type == "Mark"),
    n_recaps    = sum(event_type == "Recapture"),
    n_tags      = n_distinct(pit_tag),
    n_mortality = sum(flag_mortality),
    n_yearling  = sum(flag_yearling),
    n_concern   = sum(flag_concern),
    .groups     = "drop"
  ) |>
  arrange(
    species_run, rear, release_year, method_group, release_site, season
  )

# add the number of candidate records contributed by each source cache
input_manifest = input_manifest |>
  left_join(
    candidate_marks |>
      count(mrr_project, cache_year, name = "n_candidates"),
    by = c("project_code" = "mrr_project", "file_year" = "cache_year")
  ) |>
  mutate(n_candidates = replace_na(n_candidates, 0L))

#------------------------
# SAVE PREPARED MARK DATA

# record how mark data was prepared
prep_manifest = list(
  prep_id            = prep_id,
  created_at         = Sys.time(),
  script             = "02_prep_mark_data.R",
  project_codes      = project_codes,
  file_years         = file_years,
  pdv_option         = pdv_option,
  mortality_codes    = mortality_codes,
  yearling_codes     = yearling_codes,
  concern_codes      = concern_codes,
  species_run        = sort(unique(candidate_marks$species_run)),
  event_years        = sort(unique(candidate_marks$event_year[!is.na(candidate_marks$event_year)])),
  release_years      = sort(unique(candidate_marks$release_year[!is.na(candidate_marks$release_year)])),
  n_source_events    = sum(input_manifest$n_events),
  n_candidate_events = nrow(candidate_marks),
  source_files       = input_manifest,
  session_info       = sessionInfo()
)

# create a new output folder without overwriting a previous run
dir.create(output_dir, recursive = TRUE)

# save input & prep manifests and prepared mark data
write_csv(input_manifest, file.path(output_dir, "input_manifest.csv"))
saveRDS(prep_manifest, file.path(output_dir, "prep_manifest.rds"))
write_csv(mark_availability, file.path(output_dir, "mark_availability.csv"))
saveRDS(candidate_marks, file.path(output_dir, "candidate_marks.rds"))
message("Prepared data saved to: ", output_dir)

### END SCRIPT