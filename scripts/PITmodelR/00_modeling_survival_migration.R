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
vignette("modeling-survival-migration", package = "PITmodelR")

# load needed libraries
library(PITmodelR)
library(tidyverse)

#-----------------------------
# 1) retrieve and inspect data

# pull site observations (auto-paginates by default)
obs = get_site_observations(site_code = "LGR", year = 2023)
