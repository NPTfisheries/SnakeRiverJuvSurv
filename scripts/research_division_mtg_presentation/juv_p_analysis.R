# -----------------------
# Author(s): Mike Ackerman
# Purpose: Create tag groups and fit CJS models to evaluate juvenile detection probs for key sites
#   for Cameron & Savion's IPTDS presentation
# 
# Created Date: January 21, 2026
#   Last Modified: January 22, 2026
#
# Notes: 

#------
# setup

# clear environment
rm(list = ls())

# load libraries
library(tidyverse)
library(magrittr)
library(marked)

#-----------------
# compile mrr data

# read needed mrr_files into a list
mrr_ls = list.files(path = "./data/mrr_ptagis_retrievals/",
                    pattern = "(CDR|OGR).*\\.rda$",            # CDR = JOHTRP & SECTRP; OGR = LOSTIR
                    full.names = T,
                    recursive = T) %>%
  set_names(basename(.)) %>%
  map(~{
    e = new.env(parent = emptyenv())
    nm = load(.x, envir = e)
    e[[nm[[1]]]]
  })

# flatten and join sessions and events
safe_tbl <- function(comp) {
  purrr::possibly(~ dplyr::as_tibble(.x), otherwise = tibble())(comp)
}

mrr_df = imap_dfr(mrr_ls, 
                  ~ safe_tbl(purrr::pluck(.x, "events")) %>% 
                    mutate(rda = .y)) %>%
  left_join(
    imap_dfr(mrr_ls, 
             ~ safe_tbl(purrr::pluck(.x, "sessions")) %>% 
               mutate(rda = .y)),
    by = "file_name")

#-------------------
# ZEN & SFG Analysis

# create secesh rst marks
sct_marks = mrr_df %>%
  filter(
    release_site == "SECTRP",
    species_run_rear_type == "12W",
    pit_tag != "..........",
    migration_year >= 2018,
    month(release_date) >= 7,
    !grepl("Recapture", event_type),
    !grepl("Y", conditional_comments)
  )

# retrieve & save complete tag histories
#sct_cths = get_batch_tag_histories(tag_codes = sct_marks$pit_tag)
#save(sct_cths, file = "./data/cths/SECTRP_12W_cths.rda")

# load sct_cths
load("./data/cths/SECTRP_12W_cths.rda")

# remove a select few individuals where the mark event is not SECTRP
sct_cths %<>%
  group_by(tag_code) %>%
  filter(
    !any(event_type == "Mark" & site_code != "SECTRP")
  ) %>%
  ungroup()

# define model occasions
locs = list(
  SECTRP = c("SECTRP", "SECESR"),
  ZEN    = "ZEN",
  SFG    = "SFG",
  LGR    = c("GRJ", "GRS"),
  Down   = c(
    "GOJ",                     # Little Goose Dam
    "LMJ",                     # Lower Monumental Dam
    "ICH",                     # Ice Harbor Dam
    "CRESIS",                  # Crescent Island
    "MCJ",                     # McNary Dam
    "JDJ",                     # John Day Dam
    "LMILIS",                  # Little Miller Island
    "B2J", "BCC",              # Bonneville Dam
    "PD7", "PD8", "PDW", "TWX" # Lower Columbia & Estuary
  )
)

# build capture histories
sct_ch = build_mark_histories(
  tag_history = sct_cths,
  locs_def = locs,
  site_col = "site_code",
  tag_col = "tag_code",
  time_col = "event_date",
  enforce_order = TRUE,
  keep_unknown = FALSE
) 

# frequency of capture histories
sct_ch$ch_freq

# add migration year to sct_ch$data
mod_df = sct_ch$ch_data %>%
  left_join(sct_marks %>%
              select(pit_tag, migration_year) ,
            by = c("tag_code" = "pit_tag"))

# create model data
mod_dat = data.frame(
  ch = mod_df$ch,
  freq = rep(1, dim(mod_df)[1]),
  my = as.factor(mod_df$migration_year)
)

# evaluate models  
fitCJSmodels = function(){
  
  # apparent survivals
  Phi.time    <- list(formula = ~ time)
  Phi.my      <- list(formula = ~ my + time)
  Phi.my.time <- list(formula = ~ my * time)
  
  # detection probability
  p.time      <- list(formula = ~ time)
  p.my        <- list(formula = ~ my + time)
  p.my.time   <- list(formula = ~ my * time)
  
  # construct all combinations and put into one model table
  cml     <- create.model.list(c("Phi","p")) # makes all possible combinations of those parameter formulas
  results <- crm.wrapper(cml, 
                         data = sct_proc, 
                         ddl = sct_ddl,
                         external = FALSE, 
                         accumulate = FALSE, 
                         hessian = TRUE)
  return(results)
}
sct_proc   = process.data(mod_dat, begin.time = 1, model = "CJS", groups = "my")
sct_ddl    = make.design.data(sct_proc)
sct_models = fitCJSmodels() 

sct_models
best_mod = sct_models[[5]] # Phi(~my * time)p(~my * time)

best_res = best_mod$results$reals %>%
  bind_rows(.id = "param")

# ZEN p
zen_p = best_res %>%
  filter(param == "p",
         time == 2) %>%
  ggplot(aes(x = fct_rev(my), y = estimate)) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.4) +
  geom_point() +
  coord_flip() +
  labs(title = "Juvenile Wild Chinook Released After July 1",
       x = "Migration Year",
       y = "P(Detection ZEN)") +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0))
zen_p
ggsave(zen_p, filename = "./output/figures/research_division_mtg_presentation/ZEN_p.png")

# SFG p
sfg_p = best_res %>%
  filter(param == "p",
         time == 3) %>%
  ggplot(aes(x = fct_rev(my), y = estimate)) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.4) +
  geom_point() +
  coord_flip() +
  labs(title = "Juvenile Wild Chinook Released After July 1",
       x = "Migration Year",
       y = "P(Detection SFG)") +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0))
sfg_p
ggsave(sfg_p, filename = "./output/figures/research_division_mtg_presentation/SFG_p.png")

# SFSR Phi
sfg_phi = best_res %>%
  filter(param == "Phi",
         !time == 4,
         occ == 2
         ) %>%
  mutate(reach = case_when(
    occ == 1 ~ "SECTRP to ZEN",
    occ == 2 ~ "ZEN to SFG",
    occ == 3 ~ "SFG to LGR"
  ),
  reach = factor(reach, levels = c("SECTRP to ZEN", "ZEN to SFG", "SFG to LGR"))) %>%
  # calculate cumulative survival for each reach
  group_by(my) %>%
  mutate(cum_phi = cumprod(estimate)) %>%
  ungroup() %>%
  # start plot
  ggplot(aes(x = my, y = estimate)) +
  geom_errorbar(
    aes(ymin = lcl, ymax = ucl),
    width = 0.2,
    linewidth = 0.6
  ) +
  geom_point(size = 2) +
  #facet_wrap(~ reach, ncol = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = scales::breaks_pretty(5)) +
  labs(
    title = "ZEN to SFG, Wild Juvenile Chinook",
    x = "Migration Year",
    y = "P(Survival)"
  ) +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0))
sfg_phi  
ggsave(sfg_phi, filename = "./output/figures/research_division_mtg_presentation/SFG_phi.png")

#-------------------
# WR2 & WR1 Analysis

# create secesh rst marks
lost_marks = mrr_df %>%
  filter(
    release_site == "LOSTIR",
    capture_method == "SCREWT",
    species_run_rear_type == "11W",
    pit_tag != "..........",
    #migration_year >= 2018,
    month(release_date) >= 7,
    !grepl("Recapture", event_type),
    !grepl("Y", conditional_comments)
  )

# retrieve & save complete tag histories
lost_cths = get_batch_tag_histories(tag_codes = lost_marks$pit_tag)
save(lost_cths, file = "./data/cths/LOSTIR_11W_cths.rda")

# load sct_cths
load("./data/cths/LOSTIR_11W_cths.rda")

# remove a select few individuals where the mark event is not SECTRP
# lost_cths %<>%
#   group_by(tag_code) %>%
#   filter(
#     !any(event_type == "Mark" & site_code != "SECTRP")
#   ) %>%
#   ungroup()

# define model occasions
locs = list(
  LOSTIR = "LOSTIR",
  WR2    = "WR2",
  WR1    = "WR1",
  LGR    = c("GRJ", "GRS"),
  Down   = c(
    "GOJ",                     # Little Goose Dam
    "LMJ",                     # Lower Monumental Dam
    "ICH",                     # Ice Harbor Dam
    "CRESIS",                  # Crescent Island
    "MCJ",                     # McNary Dam
    "JDJ",                     # John Day Dam
    "LMILIS",                  # Little Miller Island
    "B2J", "BCC",              # Bonneville Dam
    "PD7", "PD8", "PDW", "TWX" # Lower Columbia & Estuary
  )
)

# build capture histories
lost_ch = build_mark_histories(
  tag_history = lost_cths,
  locs_def = locs,
  site_col = "site_code",
  tag_col = "tag_code",
  time_col = "event_date",
  enforce_order = TRUE,
  keep_unknown = FALSE
) 

# frequency of capture histories
lost_ch$ch_freq

# add migration year to sct_ch$data
mod_df = lost_ch$ch_data %>%
  left_join(lost_marks %>%
              select(pit_tag, migration_year) ,
            by = c("tag_code" = "pit_tag")) %>%
  # let's ignore migration year 2026 for now, probably not complete outmigration
  filter(!migration_year == 2026)

# create model data
mod_dat = data.frame(
  ch = mod_df$ch,
  freq = rep(1, dim(mod_df)[1]),
  my = as.factor(mod_df$migration_year)
)

# evaluate models  
fitCJSmodels = function(){
  
  # apparent survivals
  Phi.time    <- list(formula = ~ time)
  Phi.my      <- list(formula = ~ my + time)
  Phi.my.time <- list(formula = ~ my * time)
  
  # detection probability
  p.time      <- list(formula = ~ time)
  p.my        <- list(formula = ~ my + time)
  p.my.time   <- list(formula = ~ my * time)
  
  # construct all combinations and put into one model table
  cml     <- create.model.list(c("Phi","p")) # makes all possible combinations of those parameter formulas
  results <- crm.wrapper(cml, 
                         data = lost_proc, 
                         ddl = lost_ddl,
                         external = FALSE, 
                         accumulate = FALSE, 
                         hessian = TRUE)
  return(results)
}
lost_proc   = process.data(mod_dat, begin.time = 1, model = "CJS", groups = "my")
lost_ddl    = make.design.data(lost_proc)
lost_models = fitCJSmodels() 

lost_models
best_mod = lost_models[[2]] # Phi(~my * time)p(~my * time)

best_res = best_mod$results$reals %>%
  bind_rows(.id = "param")

# WR2 p
wr2_p = best_res %>%
  filter(param == "p",
         time == 2,
         !my %in% c(2018, 2024)) %>%
  ggplot(aes(x = fct_rev(my), y = estimate)) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.4) +
  geom_point() +
  coord_flip() +
  labs(title = "Juvenile Wild Chinook Released After July 1",
       x = "Migration Year",
       y = "P(Detection WR2)") +
  theme_bw() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0))
wr2_p
ggsave(wr2_p, filename = "./output/figures/research_division_mtg_presentation/WR2_p.png")

# WR2 p
# wr1_p = best_res %>%
#   filter(param == "p",
#          time == 3) %>%
#   ggplot(aes(x = fct_rev(my), y = estimate)) +
#   geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.4) +
#   geom_point() +
#   coord_flip() +
#   labs(title = "Juvenile Wild Chinook Released After July 1",
#        x = "Migration Year",
#        y = "P(Detection WR1)") +
#   theme_bw() +
#   theme(axis.title = element_text(size = 12, face = "bold"),
#         axis.text = element_text(size = 10),
#         plot.title = element_text(size = 12, face = "bold", hjust = 0))
# wr1_p
# ggsave(wr1_p, filename = "./output/figures/research_division_mtg_presentation/WR1_p.png")

# LOSTIR Phi
# lost_phi = best_res %>%
#   filter(param == "Phi",
#          !time == 4) %>%
#   mutate(reach = case_when(
#     occ == 1 ~ "LOSTIR to WR2",
#     occ == 2 ~ "WR2 to WR1",
#     occ == 3 ~ "WR1 to LGR"
#   ),
#   reach = factor(reach, levels = c("LOSTIR to WR2", "WR2 to WR1", "WR1 to LGR"))) %>%
#   # calculate cumulative survival for each reach
#   group_by(my) %>%
#   mutate(cum_phi = cumprod(estimate)) %>%
#   ungroup() %>%
#   # start plot
#   ggplot(aes(x = my, y = estimate)) +
#   geom_errorbar(
#     aes(ymin = lcl, ymax = ucl),
#     width = 0.2,
#     linewidth = 0.6
#   ) +
#   geom_point(size = 2) +
#   facet_wrap(~ reach, ncol = 1) +
#   coord_cartesian(ylim = c(0, 1)) +
#   scale_y_continuous(breaks = scales::breaks_pretty(5)) +
#   labs(
#     x = "Migration Year",
#     y = "P(Survival)"
#   ) +
#   theme_bw() +
#   theme(axis.title = element_text(size = 12, face = "bold"),
#         axis.text = element_text(size = 10),
#         plot.title = element_text(size = 12, face = "bold", hjust = 0))
# lost_phi

### END SCRIPT