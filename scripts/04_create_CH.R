# Purpose: Summarize mark events and load observation file data from PTAGIS
# Author: Ryan N. Kinzer
# Created Date: April 16, 2025


# load libraries
library(PITmodelR)
library(tidyverse)

# Get Data ---- 

tag_yr = c(2024,2025)
mig_yr = 2025

mrr_dat_ls <- readRDS(paste0('./output/mrr_files/npt_mrr_files_',gsub(', ',"_",toString(tag_yr)),'.rds'))


parr_marks <- parr_tally %>%
  filter(event_type == 'Mark') %>%
  mutate(release_year = year(release_ym)) %>%
  group_by(release_year, release_site, species_run_rear_type, capture_method) %>%
  summarise(n_tags = sum(n_tags))

ggplot(parr_marks) +
  geom_col(aes(x = release_site, y = n_tags, fill = capture_method)) +
  facet_wrap(~species_run_rear_type, scales = 'free') +
  coord_flip()

# calculate survival

kfh <- hat_events %>%
  filter(mrr_project == 'SCS',
         user_ext == 'KFH',
         year(release_ym) == 2025)

locs <- list(
  DWOR = "DWOR", # site_code is showing the event_site from the MRR file which is the mark site
  #LAP    = "LAP",
  CWR    = "CWR",
  LGR    = c("GRJ", "GRS"), # Lower Granite Juvenile Bypass (J) and Spillway (S)
  # upstream -> downstream
  Down   = c(
    "GOJ",                     # Little Goose Dam
    "LMJ",                     # Lower Monumental Dam
    "ICH",                     # Ice Harbor Dam
    "MCJ",                     # McNary Dam
    "JDJ",                     # John Day Dam
    "B2J", "BCC",              # Bonneville Dam
    "PD5", "PD6", "PD7", "PD8", "PDO", "PDW", "TWX" # Lower Columbia & Estuary
  )
)


api_key <- "1AA5CF55-C98E-4001-96E4-1E96CEE1E806"
Sys.setenv(PTAGIS_API_KEY = api_key)


tmp <- sample_n(kfh %>% filter(pit_tag != ".........."), 1000)
write_delim(tibble(tmp$pit_tag), file = './output/mrr_files/test_tags.txt')
