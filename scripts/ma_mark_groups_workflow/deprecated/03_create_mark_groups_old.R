# CREATE MARK GROUPS
#
# Author(s): Mike Ackerman, Michelle Briggs, and Ryan Kinzer
#
# Date Created: November 20, 2026
#   Last Modified: November 23, 2026
#
# Purpose: Select juvenile release groups from prepared PTAGIS MRR data,
#   apply exclusions and QC, and retain mark groups meeting a minimum sample size.
#
# NOTES:
# - Both Mark and Recapture events are eligible release events.
# - Retain the earliest eligible event per PIT tag within each mark group.
# - A PIT tag may occur in multiple groups e.g., initial marks and downstream recaps.
# - Exclusion flags and initial seasons, method_groups, etc. are defined in 02_prep_mark_data.R.

# clear environment
rm(list = ls())

# load necessary packages
library(tidyverse)
library(here)

#------------------------
# LOAD PREPARED MARK DATA
prep_id  = "prep_20261123_1259"
prep_dir = here("data", "processed", "mark_prep", prep_id)

candidate_marks = readRDS(file.path(prep_dir, "candidate_marks.rds"))

#-------------------------
# SELECT INITIAL MARK DATA

# NULL retains all available values
species_keep       = "SpSu_Chinook"
rear_keep          = "W"
release_years_keep = 2021L:2025L
methods_keep       = c("ROVING", "SCREWT")

available_marks = candidate_marks |>
  filter(
    is.null(species_keep)       | species_run %in% species_keep,
    is.null(rear_keep)          | rear %in% rear_keep,
    is.null(release_years_keep) | release_year %in% release_years_keep,
    is.null(methods_keep)       | method_group %in% methods_keep
  )

#----------------------------
# EXPLORE AVAILABLE MARK DATA
site_inventory = available_marks |>
  group_by(method_group, release_site, release_year, season) |>
  summarise(
    n_events = n(),
    n_tags = n_distinct(pit_tag),
    .groups = "drop"
  ) |>
  arrange(method_group, release_site, release_year, season)

#----------------------
# SELECT RELEASE GROUPS

# select release sites separately for each capture-method group
# NULL retains all sites for that method
sites_keep = list(
  SCREWT = c(
    "JOHTRP", "SECTRP", "SFSRKT"
  ),
  ROVING = c(
    # Dry Clearwater
    "AMERR", "CROOKR", "NEWSOC", "REDR",
    # Grande Ronde
    "CATCNF", "CATCSF", "CATHEC", "GRAND2", "HURRIC", "LOSTIR", "MINAMR", "WALLOR",
    # Imnaha
    "IMNAHR",
    # Middle Fork Salmon
    "BEARVC", "BIG2C", "CAPEHC", "CHAMBC", "ELKC", "LOONC", "MARSHC", "VALEYC",
    # South Fork Salmon
    "LAKEC", "SUMITC", "GROUSC", "SECESR", "KNOXB", "SALRSF", "JOHNSC", "BURNLC", "SAEFSF",
    # Wet Clearwater
    "LOLOC"
  )
)

# select seasons and optionally combine summer/fall
seasons_keep        = c("summer", "fall")  # NULL retains all
combine_summer_fall = FALSE

# double-check whether selected sites exist in the available_marks
if (any(lengths(sites_keep) > 0L)) {
  missing_sites = imap_dfr(
    sites_keep[lengths(sites_keep) > 0L],
    ~ tibble(method_group = .y, release_site = .x)
  ) |>
    anti_join(
      distinct(site_inventory, method_group, release_site),
      by = c("method_group", "release_site")
    )
  
  if (nrow(missing_sites) > 0) {
    warning("Selected release sites not found in available data:")
    print(missing_sites, n = Inf)
  }
}

#-----------------------
# APPLY EXCLUSIONS AND QC

# apply user-selected exclusion rules
exclude_mortality = TRUE
exclude_yearlings = TRUE  # SpSu_Chinook only, for now; see 02_prep_mark_data.R
exclude_concerns  = TRUE

# minimum number of unique PIT tags per mark group
min_tags = 100L

# variables defining each mark group
group_vars = c(
  "species_run", "rear", "release_year",
  "method_group", "release_site", "season_group"
)

# select release sites/seasons and flag records for exclusion
group_candidates = available_marks |>
  filter(
    # restrict sites separately by capture method
    map2_lgl(
      method_group, release_site,
      \(method, site) is.null(sites_keep[[method]]) ||
        site %in% sites_keep[[method]]
    ),
    is.null(seasons_keep) | season %in% seasons_keep
  ) |>
  mutate(
    # optionally pool summer and fall
    season_group = if_else(
      combine_summer_fall & season %in% c("summer", "fall"),
      "summer_fall",
      season
    ),
    # exclusion flags
    qc_mortality = exclude_mortality & coalesce(flag_mortality, FALSE),
    qc_yearling  = exclude_yearlings & coalesce(flag_yearling, FALSE),
    qc_concern   = exclude_concerns  & coalesce(flag_concern, FALSE),
    qc_missing = is.na(event_date) | is.na(release_date) |
      is.na(release_year) |
      is.na(release_site) |
      str_trim(coalesce(release_site, "")) == "" |
      is.na(season_group) |
      season_group %in% c("missing_date", "outside_season"),
    qc_exclude = qc_mortality | qc_yearling | qc_concern | qc_missing
  )

#----------------------
# SUMMARIZE MARK GROUPS

# summarize exclusions, repeated tags, and final group sizes
group_summary = group_candidates |>
  group_by(across(all_of(group_vars))) |>
  summarise(
    n_start     = n(),
    n_mortality = sum(qc_mortality),
    n_yearling  = sum(qc_yearling),
    n_concern   = sum(qc_concern),
    n_missing   = sum(qc_missing),
    n_excluded  = sum(qc_exclude),
    n_eligible  = sum(!qc_exclude),
    n_tags      = n_distinct(pit_tag[!qc_exclude]),
    .groups     = "drop"
  ) |>
  mutate(
    n_repeats_removed = n_eligible - n_tags,
    keep_group        = n_tags >= min_tags
  ) |>
  select(everything(), n_tags, keep_group) |>
  arrange(across(all_of(group_vars)))

#-------------------
# RETAIN FINAL MARKS

# keep first eligible release event per PIT tag within each group, then retain groups meeting min_tags
final_marks = group_candidates |>
  filter(!qc_exclude) |>
  arrange(event_date, release_date, file_name, sequence_number) |>
  distinct(across(all_of(c(group_vars, "pit_tag"))), .keep_all = TRUE) |>
  semi_join(
    group_summary |> filter(keep_group),
    by = group_vars
  )

### END SCRIPT