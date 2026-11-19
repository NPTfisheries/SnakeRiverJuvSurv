# CREATE MARK/TAG GROUPS
# 
# Author(s): Mike Ackerman, Michelle Briggs, and Ryan Kinzer
# 
# Date Created: November 19, 2026
#   Last Modified: 
#
# Purpose: Compile cached PTAGIS MRR data, explore variables, and prepare for
#   creating mark groups
#
# NOTES: Initial development focuses on spring/summer Chinook salmon and steelhead
#   captured in 2021-2025.

# clear environment
rm(list = ls())

# load necessary packages
library(tidyverse)
library(here)

#--------------------------------------------
# SELECT PROJECT/YEAR COMBINATIONS TO COMPILE
project_codes_of_interest  = c("CDR", "GAA", "PJC", "SCS", "BDA", "MAP")
calendar_years_of_interest = 2021L:2025L

project_years = crossing(
  project_code = project_codes_of_interest,
  year         = calendar_years_of_interest
)

#------------------------
# COMPILE CACHED MRR DATA
cache_dir = here("data", "cache", "ptagis", "mrr")

# identify cache files corresponding to selected projects/years
cache_files = project_years |>
  mutate(cache_path = file.path(
    cache_dir, project_code,
    paste0("mrr_", project_code, "_", year, ".rds")
  ))

# check for missing caches
missing_files = cache_files %>%
  filter(!file.exists(cache_path))

if (nrow(missing_files) > 0) {
  stop("Missing MRR caches:\n", paste(missing_files$cache_path, collapse = "\n"), "\nRetrieve additional project-years in 01_retrieve_mrr.R or revise project_years.")
}

# read cached RDS objects, retaining project/year provenance
mrr_ls = purrr::pmap(
  cache_files,
  function(project_code, year, cache_path) {
    
    x = readRDS(cache_path)
    
    if (!all(c("sessions", "events", "pdv_map") %in% names(x))) {
      stop("Unexpected cache structure: ", cache_path)
    }
    x
  }
)

# compile tables, preserving their separate structures
mrr = list(
  sessions = map_dfr(mrr_ls, "sessions"),
  events   = map_dfr(mrr_ls, "events"),
  pdv_map  = map_dfr(mrr_ls, "pdv_map")
)

# summarize data loaded from each cache
# cache_summary = cache_files |>
#   select(project_code, year) |>
#   mutate(
#     n_sessions = map_int(mrr_ls, ~ nrow(.x$sessions)),
#     n_events   = map_int(mrr_ls, ~ nrow(.x$events)),
#     retrieved_dttm = map_chr(mrr_ls, ~ purrr::pluck(.x, "retrieval", "retrieved_dttm", .default = NA_character_))
#   )
# print(cache_summary, n = Inf)

#-------------------
# PREPARE EVENT DATA

# confirm file_name is unique within sessions
if(anyDuplicated(mrr$sessions$file_name)) {
  stop("Duplicate file_name values found in MRR sessions.")
}

# join session information to individual events
mrr_events = mrr$events |>
  left_join(
    mrr$sessions,
    by = "file_name"
  )

# confirm join didn't change the number of events
stopifnot(nrow(mrr_events) == nrow(mrr$events))

