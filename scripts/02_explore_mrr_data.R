# -----------------------
# Author(s): Mike Ackerman
# Purpose: Compile MRR data from 01_retrieve_mrr_data.R, validate, explore, and prepare
#   for creating mark groups.
# 
# Created Date: November 18, 2025
#   Last Modified: November 20, 2025
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

# get paths to all mrr_ls .rda files
mrr_ls_files = list.files(
  path = "./output/mrr_ptagis_retrievals/",
  pattern = "\\.rda$",
  full.names = T,
  recursive = T
)

read_rda_object <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  obj_name <- ls(e)[1]   # assume exactly one object saved
  e[[obj_name]]
}

all_mrr_ls <- purrr::set_names(
  purrr::map(mrr_ls_files, read_rda_object),
  basename(mrr_ls_files)   # name each element by the .rda file for provenance
)

# helper: safely pull a component that may/may not exist
get_comp <- function(x, nm) {
  if (is.list(x) && !is.null(x[[nm]]) && is.data.frame(x[[nm]])) {
    tibble::as_tibble(x[[nm]])
  } else {
    tibble::tibble()
  }
}

sessions_df <- purrr::imap_dfr(all_mrr_ls, ~ dplyr::mutate(get_comp(.x, "sessions"), rda = .y))
events_df   <- purrr::imap_dfr(all_mrr_ls, ~ dplyr::mutate(get_comp(.x, "events"),   rda = .y))
pdv_map_df  <- purrr::imap_dfr(all_mrr_ls, ~ dplyr::mutate(get_comp(.x, "pdv_map"),  rda = .y))

# test a join
events_plus_session = events_df %>%
  left_join(
    sessions_df,
    by = "file_name"
  )

# which release sites are in each mrr project?
table(events_plus_session$release_site, events_plus_session$mrr_project)

