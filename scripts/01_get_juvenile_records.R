# Purpose: Retrieve juvenile MRR records from uploaded PTAGIS files and write
# responses from get_batch_file_data() to file.
# Author: Ryan Kinzer
# Created: 16 March 2026
# Modified: 


# load libraries
#remotes::install_github("ryankinzer/PITmodelR", ref = "main", build_vignettes = T, force = T)
library(PITmodelR)
library(tidyverse)

# parameters

tag_yr <- c(2024,2025) # could be multiple years c(2025, 2026)

proj_codes <- c('BDA', # Bill Arnsberg
                'SCS', # Sherman Sprague
                'NPC', # NPT Supplementation Evaluations
                'PJC', # Peter Cleary - Grande Ronde Supplementation/Lostine
                'IMN', # Imnahna Trap
                'CDR', # Craig Rabe
                'GAA') # Gordan Axtel - Idaho Parr Tagging


proj_years <- expand.grid('mrr_project' = proj_codes,
                          'year' = tag_yr,
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

saveRDS(mrr_dat_ls, file = paste0('./output/mrr_files/npt_mrr_files_',gsub(', ',"_",toString(tag_yr)),'.rds'))

View(mrr_dat_ls$sessions)
View(mrr_dat_ls$events)
View(mrr_dat_ls$pdv_map)
