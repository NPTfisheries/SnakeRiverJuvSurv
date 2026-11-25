# CREATE MARK GROUPS
#
# Author(s): Mike Ackerman, Michelle Briggs, Ryan Kinzer
#
# Date Created: November 24, 2026
#   Last Modified:
#
# Purpose: Use a YAML configuration with prepared MRR Mark/Recpature events.
#   Explore available sample sizes, candidate groups, apply exclusions, and 
#   duplicate tag rules, and review resulting mark groups before exporting.
#
# Notes:
# - All decisions (e.g., species, methods, seasons, exclusions) are in YAML.
# - Edit the YAML, not the processing functions, to change group definitions.
# - Mark and Recapture events are both eligible
# - Multiple release groups may legitimately contain the same PIT tag e.g., an electrofished juvenile recaptured in a screw trap.
# - Configuration lives in config/mark_groups/.

# clear environment
rm(list = ls())

# load necessary packages
library(tidyverse)
library(here)
library(yaml)

#------------------------
# LOAD CONFIGURATION FILE
config_file = here("config", "mark_groups", "chnk_parr_michelle.yml")
config      = yaml::read_yaml(config_file)

source(here("R", "validate_mark_config.R"))

# validate the configuration file
lookups = validate_mark_config(config)

#-------------------
# CREATE MARK GROUPS
source(here("R", "create_mark_groups.R"))

# load mark/recapture events for prep_id
mark_events = readRDS(here("data", "processed", "prepped_marks", config$prep_id, "mark_events.rds"))

results = create_mark_groups(mark_events, config, lookups)

message("Available mark groups: ",            nrow(results$mark_group_summary), "\n",
        "Selected mark groups: ",             sum(results$mark_group_summary$selected), "\n",
        "Selected groups meeting min_tags: ", sum(results$mark_group_summary$keep_group, na.rm = TRUE), "\n",
        "Final marks in qualifying groups: ", nrow(results$final_marks))

head(results$mark_group_summary)   # available and selected groups, with QC summary
head(results$final_marks)          # one eligible release record per PIT/group

#-----------------------------------------------
# EXPORT (ONLY WHEN mode: export IN CONFIG FILE)
source(here("R", "export_mark_groups.R"))
export_dir = export_mark_groups(results, config, config_file)

### END SCRIPT
