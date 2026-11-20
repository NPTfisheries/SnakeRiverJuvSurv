# CREATE MARK GROUPS
#
# Author(s): Mike Ackerman, Michelle Briggs, and Ryan Kinzer
#
# Date Created: November 20, 2026
#   Last Modified: 
#
# Purpose: Select candidate release groups, apply exclusions and minimum sample sizes, 
#   and QC prior to exporting mark groups.
#
# NOTES: Mark and Recapture events are both eligible to include in mark groups.

rm(list = ls())

library(tidyverse)
library(here)

#------------------------
# SELECT PREPARED MARK DATA

# identify the preparation run to use
prep_id = "prep_20261120_1616"
mark_dir = here("data", "processed", "mark_prep", prep_id)

candidate_marks = readRDS(file.path(mark_dir, "candidate_marks.rds"))
prep_manifest   = readRDS(file.path(mark_dir, "prep_manifest.rds"))

#-------------------------
# DEFINE MARK GROUPS

# select species, rear, years, methods, sites, and seasons
# NULL includes all available values for that variable
species_keep       = "SpSu_Chinook"
rear_keep          = "W"
release_years_keep = 2021L:2025L
methods_keep       = c("ROVING", "SCREWT")
sites_keep         = NULL
seasons_keep       = c("summer", "fall")

# exclusion rules
exclude_mortality = TRUE
exclude_yearlings = TRUE
exclude_concerns  = TRUE

# minimum number of unique PIT tags per group
min_tags = 100L

# variables defining individual release groups
group_vars = c(
  "species_run", "rear", "release_year",
  "method_group", "release_site", "season"
)

#----------------------------
# FILTER CANDIDATE MARK DATA

# select the biological groups and release events of interest
eligible_marks = candidate_marks |>
  filter(
    is.null(species_keep)         | species_run %in% species_keep,
    is.null(rear_keep)            | rear %in% rear_keep,
    is.null(release_years_keep)   | release_year %in% release_years_keep,
    is.null(methods_keep)         | method_group %in% methods_keep,
    is.null(sites_keep)           | release_site %in% sites_keep,
    is.null(seasons_keep)         | season %in% seasons_keep,
    # apply user-selected biological exclusions
    !(exclude_mortality & flag_mortality),
    !(exclude_yearlings & flag_yearling),
    !(exclude_concerns  & flag_concern),
    # require valid release-group information
    !is.na(release_year),
    !is.na(release_site),
    release_site != "",
    !season %in% c("missing_date", "outside_season")
  )

#-----------------------
# SUMMARIZE MARK GROUPS

# summarize qualifying events by release group
group_summary = eligible_marks |>
  group_by(across(all_of(group_vars))) |>
  summarise(
    n_events = n(),
    n_marks = sum(event_type == "Mark"),
    n_recaps = sum(event_type == "Recapture"),
    n_tags = n_distinct(pit_tag),
    .groups = "drop"
  ) |>
  mutate(
    keep_group = n_tags >= min_tags
  ) |>
  arrange(
    species_run, rear, release_year,
    method_group, release_site, season
  )

# retain individual release events from groups meeting min_tags
selected_marks = eligible_marks |>
  semi_join(
    group_summary |> filter(keep_group),
    by = group_vars
  )

# CONTINUE HERE

### END SCRIPT
