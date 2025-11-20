# -----------------------
# Author(s): Mike Ackerman
# Purpose: Compile mark groups
# 
# Created Date: November 18, 2025
#   Last Modified: 
#
# Notes: 

#--------
# Setup

# clear environment
rm(list = ls())

# install PITmodelR, if needed
#remotes::install_github("ryankinzer/PITmodelR", ref = "ma_develop", build_vignettes = T, force = T)

# load libraries
library(PITmodelR)
library(tidyverse)

#---------------
# Project Data

# all project codes of interest
proj_codes = c("AAB", # Alan Byrne Projects
               "CDR", # Craig Rabe Projects
               "GAA", # Gordon Axel Projects
               "IMN", # Nez Perce Tribe Imnaha Basin Project
               "JLV", # Jason Vogel Projects
               "KAA", # Kim Apperson Projects
               "MAP", # McCall Anadromous Projects
               "NPC", # NPT Supplementation Eval
               "NPM", # Idaho Natural Production Monitoring and Evaluation Project
               "PAK", # Paul Kucera Projects
               "RNK"  # Ryan Kinzer Projects
  ) %>%
  # remove project codes with no .xml files
  setdiff(c("AAB", "JLV", "RNK"))
  
# get years for each proj_codes
proj_yrs = map_dfr(proj_codes, function(cd) {
  yrs = get_project_years(cd)
  tibble(code = cd, year = yrs)  
}) %>%
  # just deal with more recent data, for now
  filter(year >= 2020)

# retrieve mrr files for each project code and year
all_files = proj_yrs %>%
  mutate(files = map2(code, year, ~ {
    out = get_mrr_files(.x, .y)
    out = out %>% select(-year)  # drop internal year column from ptagis
    out
  })) %>%
  unnest(files) %>%
  # duplicative with projectCode
  select(-code)

# check file types
table(all_files$projectCode, all_files$fileTypeExtension)
table(all_files$year, all_files$fileTypeExtension)

# let's just deal with .xml files, for now
mrr_files = all_files %>%
  filter(fileTypeExtension == ".xml")

#---------------------------------------------------------
# Safely Retrieve MRR Data for Each proj_yrs in mrr_files

# create directory to save output
out_dir = "./output/mrr_data"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

walk(seq_len(nrow(proj_yrs)), function(i) {
  
  cd = proj_yrs$code[i]
  yr = proj_yrs$year[i]
  
  message("Retrieving data for MRR project code ", cd, ", year ", yr)
  
  # filter files for each combo
  py_files = mrr_files %>%
    filter(projectCode == cd, year == yr)
  
  # tryCatch to avoid breaking on errors
  mrr_ls = tryCatch({
    get_batch_file_data(
      filenames = py_files$name,
      keep_code_cols = FALSE,
      label_conflict = "suffix",
      use_codes_on_conflict = TRUE
    )
  }, error = function(e) {
    message("ERROR: Failed for projectCode ", cd, " year", yr)
    message("Details: ", e$message)
    return(NULL) # continue to next
  })
  
  # skip saving if result is NULL
  if (is.null(mrr_ls)) return(NULL)
  
  # save the result
  sub_dir = file.path(out_dir, cd) 
  if (!dir.exists(sub_dir)) dir.create(sub_dir, recursive = TRUE)
  save(mrr_ls, file = file.path(sub_dir, paste0("mrr_", cd, "_", yr, ".rda")))
  message("Saved ", out_file)
})

#---------------------------------------------------------
# Compile All MRR Lists from Above

# path to mrr_ls files
mrr_ls_files = list.files(
  path = "./output/mrr_data",
  pattern = "\\.rda",
  full.names = T,
  recursive = T
)

all_mrr_ls = lapply(mrr_ls_files, function(f) {
  env = new.env()
  load(f, envir = env)
  obj = ls(env)[1]
  get(obj, envir = env)
})

combined_mrr = list(
  sessions = bind_rows(lapply(all_mrr_ls, `[[`, "sessions")),
  events   = bind_rows(lapply(all_mrr_ls, `[[`, "events")),
  issues   = bind_rows(lapply(all_mrr_ls, `[[`, "issues"))
)

### END SCRIPT
