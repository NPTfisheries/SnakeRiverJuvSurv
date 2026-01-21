# -----------------------
# Author(s): Mike Ackerman
# Purpose: Compile MRR data from 01_retrieve_mrr_data.R, validate, explore, and prepare
#   for creating mark groups.
# 
# Created Date: November 18, 2025
#   Last Modified: January 5, 2026
#
# Notes: 

#--------
# Setup

# clear environment
rm(list = ls())

# load libraries
library(PITmodelR)
library(tidyverse)

#---------------
# MRR Sites of Interest
mrr_sites_of_interest = tribble(
  ~location, ~capture_method, ~udl, ~release_site,
  "Lake Creek",    "SCREWT",                       "LCT", NA,                    # before 2020
  "Upper Secesh",  "SCREWT",                       "SRT", NA,                    # before 2020
  "Upper Secesh",  c("BSEINE", "PSEINE", "SHOCK"), NA,    c("LAKEC", "SECESR"),  # GAA     
  "Lower Secesh",  "SCREWT",                       "SCT", c("SECTRP", "SECESR"), # CDR
  "Lick Creek",    "SCREWT",                       NA,    NA,                    # before 2020
  "Upper SFSR",    "SCREWT",                       NA,    NA,                    # before 2020
  "Upper SFSR",    c("BSEINE", "PSEINE", "SHOCK"), NA,    c("KNOXB", "SALRSF"),  # GAA
  "Johnson Creek", "SCREWT",                       "JCT", c("JOHTRP", "JOHNSC"), # CDR
  "SFSR Krassel",  "SCREWT",                       NA,    "SFSRKT",              # MAP
  "Big Creek",     c("BSEINE", "PSEINE", "SHOCK"), NA,    "BIG2C",               # GAA
  "Big Creek",     "SCREWT",                       NA,    "BIG2CT",              # MAP
  "Newsome Creek", "SCREWT",                       NA,    NA,                    # SCS?
  "SF Clearwater", "SCREWT",                       NA,    "SFCTRP",              # NPC
  "Lolo Creek",    "SCREWT",                       NA,    "LOLTRP",              # SCS?
  "Imnaha River",  "SCREWT",                       NA,    "IMNTRP"               # IMN
)

#---------------------------------------------------------
# Compile MRR Outputs from 01_retrieve_mrr_data.R

# get paths to all mrr data retrieval .rda files
mrr_files = list.files(
  path = "./data/mrr_ptagis_retrievals/",
  pattern = "\\.rda$",
  full.names = T,
  recursive = T
)

# read all .rda files
all_mrr = mrr_files %>%
  set_names(basename(.)) %>%
  map(~{
    e = new.env(parent = emptyenv())
    nm = load(.x, envir = e)
    e[[nm[[1]]]]
  })

safe_tbl <- function(comp) {
  purrr::possibly(~ dplyr::as_tibble(.x), otherwise = tibble())(comp)
}

# compile sessions, events, and (optionally) pdv mapping
sessions_df <- imap_dfr(all_mrr, ~ safe_tbl(purrr::pluck(.x, "sessions")) %>% mutate(rda = .y))
events_df   <- imap_dfr(all_mrr, ~ safe_tbl(purrr::pluck(.x, "events"))   %>% mutate(rda = .y))
#pdv_map_df  <- imap_dfr(all_mrr, ~ safe_tbl(purrr::pluck(.x, "pdv_map"))  %>% mutate(rda = .y))

# join session info to events
events_session_df = events_df %>%
  left_join(
    sessions_df,
    by = "file_name"
  )

# which release sites are in each mrr project?
table(events_session_df$release_site, events_session_df$mrr_project)
table(events_session_df$mrr_project, events_session_df$capture_method)

### END SCRIPT