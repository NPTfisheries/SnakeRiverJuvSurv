# -----------------------
# Author(s): Mike Ackerman
# Purpose: Retrieve MRR data from PTAGIS for project codes and years of interest
#   and write responses from get_batch_file_data() to file.
# 
# Created Date: November 18, 2025
#   Last Modified: November 20, 2025
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

#---------------
# Project Data

# project codes of interest
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
               "RNK", # Ryan Kinzer Projects
               "SCS"  # Sherman Sprague Projects
  ) %>%
  # remove some project codes with no data for mrr_sites_of_interest (2020 or after)
  setdiff(c("AAB", "JLV", "KAA", "NPM", "PAK", "RNK"))
  
# get years for each proj_codes
proj_yrs = map_dfr(proj_codes, function(cd) {
  yrs = get_project_years(cd)
  tibble(code = cd, year = yrs)  
}) %>%
  # just deal with more recent data, for now
  filter(year >= 2020)

# retrieve mrr file info for each project code and year
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
out_dir = "./output/mrr_ptagis_retrievals"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# loop over each code x year and save results
walk(seq_len(nrow(proj_yrs)), function(i) {
  
  cd = proj_yrs$code[i]
  yr = proj_yrs$year[i]
  
  message("Retrieving data for MRR project code ", cd, ", year ", yr)
  
  # filter files for each project and year
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
  }, 
  # errors largely resolved by bug fix to check_pdv_label_consistency()
  error = function(e) {
    message("ERROR: Failed for projectCode ", cd, " year", yr)
    message("Details: ", e$message)
    return(NULL) # continue to next
  })
  
  # skip saving if result is NULL
  if (is.null(mrr_ls)) return(NULL)
  
  # save the result
  sub_dir = file.path(out_dir, cd) 
  if (!dir.exists(sub_dir)) dir.create(sub_dir, recursive = TRUE)
  out_file = file.path(sub_dir, paste0("mrr_", cd, "_", yr, ".rda"))
  save(mrr_ls, file = out_file)
  message("Saved ", out_file)
})

### END SCRIPT
