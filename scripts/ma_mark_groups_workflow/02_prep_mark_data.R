# PREPARE MRR MARK DATA
#
# Author(s): Mike Ackerman, Michelle Briggs, Ryan Kinzer
#
# Date Created: November 24, 2026
#   Last Modified:
#
# Purpose: Compile cached PTAGIS MRR Mark/Recapture events into one "analysis-neutral"
#   dataset used for mark-group creation.
#
# Notes:
# - 01_retrieve_mrr.R retrieves and caches; this script does not.
# - No species, rear, capture-method, site, seasonal, size, or conditional_comment exclusions occur here.
# - Invalid/missing data and repeated PIT tags are still retained for later QC.
# - cache_year is the PTAGIS file/session/submission year, not a fish's release_year.

# clear environment
rm(list = ls())

# load necessary packages
library(tidyverse)
library(here)

#-----------------------
# SELECT CACHED MRR DATA

# NULL includes all available project codes and/or cache years
project_codes = NULL         # e.g., "CDR" or c("BDA", "CDR", "GAA")
cache_years   = NULL         # e.g., 2021L:2025L

#--------------------
# PREPARE MARK EVENTS
source(here("R", "compile_mrr_mark_events.R"))
prepared = compile_mrr_mark_events(
  cache_dir     = here("data", "cache", "ptagis", "mrr"),
  project_codes = project_codes,
  cache_years   = cache_years
)

mark_events     = prepared$mark_events
source_manifest = prepared$source_manifest

#------------------
# SAVE PREPPED DATA
prep_id    = paste0("prep_", format(Sys.time(), "%Y%m%d_%H%M"))
output_dir = here("data", "processed", "prepped_marks", prep_id)

prep_manifest = list(
  prep_id             = prep_id,
  created_at          = Sys.time(),
  script              = "02_prep_mark_data.R",
  project_codes       = project_codes,
  cache_years         = cache_years,
  event_types         = c("Mark", "Recapture"),
  n_total_mrr_events  = sum(source_manifest$n_events),
  n_mark_recap_events = nrow(mark_events),
  source_files        = source_manifest,
  session_info        = sessionInfo()
)

dir.create(output_dir, recursive = TRUE)
readr::write_csv(source_manifest, file.path(output_dir, "source_manifest.csv"))
saveRDS(prep_manifest,            file.path(output_dir, "prep_manifest.rds"))
saveRDS(mark_events,              file.path(output_dir, "mark_events.rds"))

message("Prepared ", nrow(mark_events), " Mark/Recapture events from ", nrow(source_manifest), " MRR caches: ", output_dir)

### END SCRIPT
