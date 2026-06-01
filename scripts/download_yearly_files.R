# Purpose: Compile mark file data from a single year for ESA reporting purposes.
# Author: Ryan N. Kinzer
# Created Date: December 3, 2025

# load libraries
library(PITmodelR)
library(tidyverse)

# Get Data ---- 

yr = c(2024, 2025)

all_project <- get_project_codes()



# project codes of interest
proj_codes = c("CDR", # Craig Rabe Projects
               "IMN", # Nez Perce Tribe Imnaha Basin Project
               "NPC", # NPT Supplementation Eval, Hook and Line
               "PJC", # e-fishing in Grande Ronde
               "BDA", # fall chinook seining, parr seining and e-fishing
               "SCS"  # Sherman Sprague Projects - Lostine Coho, Lolo RST, some parr seiniing/e-fishing
)

proj_years <- expand.grid('mrr_project' = proj_codes,
                          'year' = yr,
                          KEEP.OUT.ATTRS = FALSE,
                          stringsAsFactors = FALSE)

# retrieve mrr file info for each project code and year

safe_get_mrr <- possibly(
  ~ get_mrr_files(code = .x, year = .y),
  otherwise = NULL
)

all_files_ls <- map2(
  proj_years$mrr_project,
  proj_years$year,
  function(x, y) {
    out <- safe_get_mrr(x, y)
    if (is.null(out)) warning("Warning: Failed for project code: ", x, " and year: ", y)
    out
  }
  ) %>%
  compact()

all_files <- all_files_ls %>%
  bind_rows() %>%
  mutate(user_ext = str_extract(name, "(?<=-)([A-Za-z0-9]{3})(?=\\.)"))

# check file types
table(all_files$projectCode, all_files$fileTypeExtension)
table(all_files$user_ext, all_files$projectCode)

# get mrr file data

mrr_dat_ls <- get_batch_file_data(all_files$name,
                                  pdvs = "attach",
                                  map_pdvs = FALSE)

View(mrr_dat_ls$sessions)
View(mrr_dat_ls$events)
View(mrr_dat_ls$pdv_map)

# create directory to save output
out_dir = paste0("./output/mrr_files")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

saveRDS(mrr_dat_ls, paste0(out_dir,"/npt_",gsub(", ","_",toString(yr)),".rds"))

