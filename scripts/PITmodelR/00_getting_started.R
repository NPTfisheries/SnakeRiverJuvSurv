# -----------------------
# Author(s): Mike Ackerman
# Purpose: Initial testing of ryankinzer/PITmodelR
# 
# Created Date: November 4, 2025
#   Last Modified:
#
# Notes: 

#--------------------------
# Install and Load Packages

# clear environment
rm(list = ls())

# install PITmodelR, if needed
remotes::install_github("ryankinzer/PITmodelR", ref = "ma_develop", build_vignettes = T, force = T)

# view vignettes
vignette("getting-started", package = "PITmodelR")

# load needed libraries
library(PITmodelR)
library(tidyverse)

#-------------
# Project Data

project_codes = get_project_codes()                 # available project codes in ptagis
project_yrs   = get_project_years(code = "CDR")     # available years for a project
files_df = get_mrr_files(code = "CDR", year = 2024) # df of mrr files for given project and year

# get data for one mrr file
mrr_file = get_file_data(files_df$name[1])

# view objects in mrr_file
mrr_file$session
mrr_file$events
mrr_file$session_pdv_fields
mrr_file$detail_pdv_fields

# flatten all of above into single df
mrr_df = flatten_mrr_file(mrr_file)

# secesh trap 2024 file names
sct_24_filenames = files_df$name %>%
  str_subset("SCT")

# download all secesh trap 2024 files
sct_24_ls = get_batch_file_data(
  sct_24_filenames,
  check_labels = "warn",       # checks if user defined fields have the same labels, can use "error" to stop on mismatches
  keep_code_cols = FALSE,      # should we keep PDV/SPDV code columns and the user defined labels
  label_conflict = "suffix",   # behavior for file label column collisions
  use_codes_on_conflict = TRUE # prefer consistent code columns across files
)

names(sct_24_ls$files)                  # list of all individual files
sct_24_ls$files$`CDR-2024-101-SCT.xml`  # look at data for a single file

# list mismatch issues
sct_24_ls$issues

# combine session and event data into a single df
sct_24_df = sct_24_ls$combined

#------------------------
# Set Up Survival Studies

# create mark group
mark_group = sct_24_df %>%
  filter(species_run_rear_type == "12W",
         migration_year == 2025,
         between(release_date, ymd(20240901), ymd(20241231)),
         release_site == "SECTRP",
         pittag != "..........",
         !grepl("RE", text_comments),
         !grepl("Y", conditional_comments))

# get tag history for a single tag
tag_history = get_tag_history(tag_code = mark_group$pittag[1])

# get tag history for a list of tags
tag_history = get_batch_tag_histories(tag_codes = mark_group$pittag)

# what kind of observations?
table(tag_history$event_type)

# locations of observations?
table(tag_history$site_code)

# set up cjs locations
locs = c("SECTRP", "ZEN", "SFG", "LGR", "Down")
locs = list(
  SECTRP = "SECTRP",
  ZEN    = "ZEN",
  SFG    = "SFG",
  LGR    = c("GRJ", "GRS"),
  Down   = c("LMN","MCN","BON", "B2J", "BCC", "GOJ", "ICH", "JDJ", "LMJ", "MCJ", "PD7", "PD8", "PDW", "TWX")
)

res = build_mark_histories(
  tag_history = tag_history,
  locs_def    = locs,
  site_col    = "site_code",
  tag_col     = "tag_code",
  time_col    = "event_time",
  enforce_order = T,
  keep_unknown = F
)

# save, view some things from res
res$dropped_summary
res$mapping

#--------------
# Fit CJS Model
library(marked)

fit = fit_marked_cjs(ch_data = res$ch_data,
                     phi_formula = ~ time,
                     p_formula   = ~ time,
                     conf_level  = 0.95)

head(fit$phi) # survival
head(fit$p)   # detection probability

# simple plots
print(fit$plots$phi)
print(fit$plots$p)

# full model summary
summary(fit$model)

# cumulative survival
fit$cum_phi
fit$covariance_mode      # "full" if covariance used, otherwise "independence_fallback"
print(fit$plots$cum_phi)

#
proc = process.data(res$ch_data, model = "CJS")
ddl  = make.design.data(proc)

Phi.time <- list(formula = ~ time)
p.time   <- list(formula = ~ time)

mod <- crm(proc, ddl, model.parameters = list(Phi = Phi.time, p = p.time))
summary(mod)

mod$results
mod$results$reals$Phi

est_preds = predict(mod) %>% 
  map(.f = as_tibble)
