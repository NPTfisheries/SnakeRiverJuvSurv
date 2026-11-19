# RETRIEVE PTAGIS MRR DATA
# 
# Author(s): Mike Ackerman, Michelle Briggs, and Ryan Kinzer
# 
# Date Created: November 18, 2026
#   Last Modified: 
#
# Purpose: Inventory, retrieve, and cache PTAGIS MRR data by project and file year
#
# NOTES: The script deliberately does NOT filter by species, event types, capture methods
#   mark/release sites, life stage, fish size, dates/seasons, etc. Those decisions
#   belong to later scripts and data preparation to create mark groups.

# clear environment
rm(list = ls())

# install PITmodelR, if needed
#remotes::install_github("ryankinzer/PITmodelR", ref = "dev", build_vignettes = F, force = T)

# load necessary packages
library(tidyverse)
library(here)
library(PITmodelR)

#---------------------------------------------
# SELECT PROJECT/YEAR COMBINATIONS TO RETRIEVE
project_codes_of_interest = c(
  "CDR",                        # Craig Rabe Projects
  "GAA",                        # Gordon Axel Projects (NOAA)
  "PJC",                        # Grande Ronde Supplementation LR M&E (NPT NE OR)
  "SCS",                        # Sherman Sprague Projects (NPT OFO)
  "BDA",                        # Bill Arnsberg Projects (NPT OFO)
  "MAP"                         # McCall Anadromous Projects
)

calendar_years_of_interest = 2025L:2026L

project_years = crossing(
  project_code  = project_codes_of_interest, 
  year          = calendar_years_of_interest
)

# optionally, create a custom tibble
# project_years = tribble(
#   ~project_code, ~year,
#   "CDR", 2025L,
# )

#--------------
# CACHE OPTIONS

# where would you like to cache the retrieved MRR data by project-year?
cache_dir = here("data", "cache", "ptagis", "mrr")

# IMPORTANT: do you want to skip caches that already exist (FALSE) or re-download and update existing caches (TRUE)?
refresh = FALSE

#--------------------------
# RETRIEVE & CACHE MRR DATA

# load cache_mrr_data() function
source(here("R", "cache_mrr_data.R"))

# retrieve each project/year and inventory
retrieval_summary = map2_dfr(
  .x = project_years$project_code,
  .y = project_years$year,
  ~ cache_mrr_data(
    project_code = .x,
    year         = .y,
    cache_dir    = cache_dir,
    refresh      = refresh
  )
)

#-------------------------
# UPDATE MRR RETRIEVAL LOG
log_file = file.path(cache_dir, "mrr_retrieval_log.csv")

# retain new info for new retrieval attempts (excludes previously cached data)
new_retrievals = retrieval_summary |>
  filter(status != "cached")

# append results to running log of MRR retrievals, if there are new retrievals
if (nrow(new_retrievals) > 0) { write_csv(new_retrievals, log_file, append = file.exists(log_file)) }

### END SCRIPT