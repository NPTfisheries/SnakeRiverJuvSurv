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

#---------------------------------------------------------
# Compile MRR Lists from 01_retrieve_mrr_data.R

# get paths to all mrr_ls .rda files
mrr_ls_files = list.files(
  path = "./output/mrr_ptagis_retrievals/",
  pattern = "\\.rda",
  full.names = T,
  recursive = T
)

all_mrr_ls = map(mrr_ls_files, function(f) {
  e = new.env()
  load(f, envir = e)
  e[[ls(e)[1]]]
})

# combine all mrr_ls files into a single list
# all_mrr_ls = lapply(mrr_ls_files, function(f) {
#   env = new.env()
#   load(f, envir = env)
#   obj = ls(env)[1]
#   get(obj, envir = env)
# })

# combine into a single list
# combined_mrr = list(
#   sessions = bind_rows(lapply(all_mrr_ls, `[[`, "sessions")),
#   events   = bind_rows(lapply(all_mrr_ls, `[[`, "events")),
#   issues   = bind_rows(lapply(all_mrr_ls, `[[`, "issues"))
# )

# let's just deal with mrr events for now
events_df = bind_rows(lapply(all_mrr_ls, `[[`, "events"))

# which release sites are in each mrr project?
table(events_df$release_site, events_df$mrrproject)

