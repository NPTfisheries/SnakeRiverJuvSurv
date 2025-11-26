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

# compile mrr_ls_files
all_mrr_ls = map(mrr_ls_files, function(f) {
  e = new.env()
  load(f, envir = e)
  e[[ls(e)[1]]]
})

# diagnosing column type conflicts
type_map <- map(all_mrr_ls, ~ .x$events) %>% 
  imap_dfr(function(df, name) {
    tibble(
      file = name,
      col  = names(df),
      class = purrr::map_chr(df, ~ paste(class(.x), collapse = ","))
    )
  })

type_conflicts <- type_map %>%
  group_by(col) %>%
  summarise(n_types = n_distinct(class),
            classes = paste(unique(class), collapse = " | "),
            .groups = "drop") %>%
  filter(n_types > 1)

# ignore sessions for now
all_events_df = map(all_mrr_ls, "events") %>% bind_rows()
all_issues_df = map(all_mrr_ls, "issues") %>% bind_rows()


# let's just deal with mrr events for now
events_df = bind_rows(lapply(all_mrr_ls, `[[`, "events"))

# which release sites are in each mrr project?
table(events_df$release_site, events_df$mrrproject)

