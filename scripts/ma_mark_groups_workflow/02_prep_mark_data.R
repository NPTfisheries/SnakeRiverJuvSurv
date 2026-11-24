# PREPARE MRR MARK DATA
#
# Author(s): Mike Ackerman, Michelle Briggs, Ryan Kinzer
#
# Date Created: November 24, 2026
#   Last Modified:
#
# Purpose: Compile cached PTAGIS MRR Mark and Recapture events into one
#   analysis-neutral dataset for mark-group creation.
#
# Notes:
# - 01_retrieve_mrr.R creates the caches; this script never retrieves PTAGIS data.
# - No species, rear, capture-method, site, season, size, mortality exclusions, etc. 
#     are applied here. Those decisions belong in the configuration YAML.
# - Repeated PIT tags and records with missing/invalid values are retained
#     for now.
# - Cache year is the PTAGIS file/session/upload year.

# clear environment
rm(list = ls())

# load necessary packages
library(tidyverse)
library(here)

#-----------------------
# SELECT CACHED MRR DATA
cache_dir = here("data", "cache", "ptagis", "mrr")

# NULL includes all available project codes and/or cache years
project_codes = NULL  # e.g., c("BDA", "CDR", "GAA")
cache_years   = NULL  # e.g., 2021L:2025L

# automatically include all matching project/year caches created by 01_retrieve_mrr.R.
cache_files = list.files(
  cache_dir,
  pattern    = "^mrr_[A-Z0-9]{3}_[0-9]{4}\\.rds$",
  recursive  = TRUE,
  full.names = TRUE
) |> sort()

if (length(cache_files) == 0L) stop("No MRR cache files found in: ", cache_dir)

cache_index = tibble(
  cache_file = normalizePath(cache_files),
  cache_project = str_match(basename(cache_files), "^mrr_([A-Z0-9]{3})_([0-9]{4})\\.rds$")[, 2],  
  cache_year = as.integer(str_match(basename(cache_files), "^mrr_([A-Z0-9]{3})_([0-9]{4})\\.rds$")[, 3])
) |>
  filter(
    is.null(project_codes) | cache_project %in% project_codes,
    is.null(cache_years)   | cache_year %in% cache_years
  )

if (nrow(cache_index) == 0L) {
  stop("No MRR caches match the selected project_codes and cache_years.")
}

# warn when a requested project or year has no selected cache.
if (!is.null(project_codes)) {
  missing_projects = setdiff(project_codes, cache_index$cache_project)
  if (length(missing_projects)) {
    warning("No matching caches for project(s): ", paste(missing_projects, collapse = ", "))
  }
}
if (!is.null(cache_years)) {
  missing_years = setdiff(cache_years, cache_index$cache_year)
  if (length(missing_years)) {
    warning("No matching caches for year(s): ", paste(missing_years, collapse = ", "))
  }
}

cache_index = cache_index |> mutate(md5 = unname(tools::md5sum(cache_file)))

message("Preparing ", nrow(cache_index), " cached project/year datasets.")

#------------------------
# COMPILE RELEASE EVENTS

# keep session attributes and original event columns, including raw PDVs. 
# join within each cache so similarly named files from different caches cannot accidentally match one another.
cache_results = map(seq_len(nrow(cache_index)), function(i) {
  path = cache_index$cache_file[i]
  mrr  = readRDS(path)
  
  if (!is.list(mrr) || !is.data.frame(mrr$sessions) || !is.data.frame(mrr$events)) {
    stop("Unexpected cache structure: ", path)
  }
  
  sessions = as_tibble(mrr$sessions)
  events   = as_tibble(mrr$events)
  
  if (!"file_name" %in% names(sessions) || !all(c("file_name", "event_type", "pit_tag", "event_date", "release_date", 
                                                  "species_run_rear_type", "capture_method", "release_site") %in% names(events))) {
    stop("Required MRR fields are missing: ", path)
  }
  
  if (anyNA(sessions$file_name) || any(!nzchar(trimws(sessions$file_name))) || anyDuplicated(sessions$file_name)) {
    stop("Missing or duplicate session file_name: ", path)
  }
  if (anyNA(events$file_name) || !all(events$file_name %in% sessions$file_name)) {
    stop("Events do not match session files: ", path)
  }
  
  # Mark and Recapture are both potential release occasions.
  releases = events |>
    filter(event_type %in% c("Mark", "Recapture")) |>
    left_join(
      sessions, 
      by = "file_name", 
      relationship = "many-to-one",
      suffix = c("", "_session")
    ) |>
    mutate(
      cache_project = cache_index$cache_project[i],
      cache_year    = cache_index$cache_year[i],
      event_year    = lubridate::year(event_date),
      release_year  = lubridate::year(release_date)
    )
  
  list(
    events            = releases,
    session_files     = sessions$file_name,
    n_sessions        = nrow(sessions),
    n_events          = nrow(events),
    n_release_events  = nrow(releases),
    n_pdv_definitions = if (is.data.frame(mrr$pdv_map)) nrow(mrr$pdv_map)
    else NA_integer_
  )
})

# guard against processing the same PTAGIS MRR file in multiple caches.
if (anyDuplicated(unlist(map(cache_results, "session_files")))) {
  stop("Duplicate MRR file_name found across caches. Check source selection.")
}

mark_events = list_rbind(map(cache_results, "events"))
if (nrow(mark_events) == 0L) stop("No Mark or Recapture events found.")

source_manifest = cache_index |>
  mutate(
    n_sessions        = map_int(cache_results, "n_sessions"),
    n_events          = map_int(cache_results, "n_events"),
    n_release_events  = map_int(cache_results, "n_release_events"),
    n_pdv_definitions = map_int(cache_results, "n_pdv_definitions")
  )

stopifnot(nrow(mark_events) == sum(source_manifest$n_release_events))

#------------------------------
# SUMMARIZE SOURCE AVAILABILITY

# inventory uses original PTAGIS codes; grouping rules are applied later.
event_inventory = mark_events |>
  count(cache_project, cache_year, species_run_rear_type,
        capture_method, release_site, release_year, event_type,
        name = "n_events") |>
  arrange(cache_project, cache_year, species_run_rear_type,
          capture_method, release_site, release_year, event_type)

#-------------------
# SAVE PREPARED DATA
prep_id    = paste0("prep_", format(Sys.time(), "%Y%m%d_%H%M"))
output_dir = here("data", "processed", "mark_prep", prep_id)

prep_manifest = list(
  prep_id         = prep_id,
  created_at      = Sys.time(),
  script          = "02_prep_mark_data.R",
  project_codes   = project_codes,
  cache_years     = cache_years,
  event_types     = c("Mark", "Recapture"),
  n_source_events = sum(source_manifest$n_events),
  n_mark_events   = nrow(mark_events),
  source_files    = source_manifest,
  session_info    = sessionInfo()
)

dir.create(output_dir, recursive = TRUE)

# write results to files
saveRDS(mark_events, file.path(output_dir, "mark_events.rds"))
write_csv(source_manifest, file.path(output_dir, "source_manifest.csv"))
write_csv(event_inventory, file.path(output_dir, "event_inventory.csv"))
saveRDS(prep_manifest, file.path(output_dir, "prep_manifest.rds"))

message("Prepared ", nrow(mark_events), " Mark/Recapture events from ", nrow(source_manifest), " caches: ", output_dir)

### END SCRIPT
