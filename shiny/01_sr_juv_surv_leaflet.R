# -----------------------
# Author: Mike Ackerman & Michelle Briggs
# Purpose: Create leaflet map to support 2025 juvenile survival work.
# 
# Created: November 5, 2026
#   Last Modified:
# 
# NOTES: 

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(here)
library(readxl)
library(PITcleanr)
library(sf)
library(leaflet)
library(htmlwidgets)

# ----------------------
# COMPILE NECESSARY DATA

# spring/summer Chinook populations
load(here("shiny", "shiny_data", "SR_pops.rda")) ; rm(fall_pop, sth_pop)

# mark release site table from multiple_release_group_simulation_estimates.Rmd
load(here("shiny", "shiny_data", "n_rel_site.rda"))

# MRR & INT sites from PTAGIS
mrr_meta = PITcleanr::queryMRRMeta()
int_meta = PITcleanr::queryInterrogationMeta()

# streams and steelhead major/minor spawning areas
load(here("shiny", "shiny_data", "steelhead_gis_data.rda"))
streams = sthd_critical %>% st_transform("EPSG:4326")

# ----------------------
# PREP DATA FOR MAPPING

# juvenile release site locations
juv_mrr_sites = mrr_meta |>
  filter(siteCode %in% unique(n_rel_site$`Release site`)) |>
  select(
    site_code = siteCode,
    name,
    rkm,
    latitude,
    longitude
  ) |>
  left_join(
    n_rel_site |>
      pivot_wider(
        id_cols = `Release site`,
        names_from = `Tag season`,
        values_from = `Number of tagged fish`
      ),
    by = c("site_code" = "Release site")
  ) |>
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )

# summarize tagged fish across seasons for site popups
site_totals = n_rel_site |>
  group_by(`Release site`) |>
  summarise(
    n_tagged = sum(`Number of tagged fish`, na.rm = TRUE),
    .groups = "drop"
  )

juv_mrr_sites = juv_mrr_sites |>
  left_join(
    site_totals,
    by = c("site_code" = "Release site")
  )

# IPTDS sites
iptds_sf = int_meta |>
  select(
    site_code = siteCode,
    name,
    rkm,
    latitude,
    longitude
  ) |>
  filter(site_code %in% c(
    "CRA", "SC1", "SC2", "SC3", "SC4",   # south fork clearwater
    "ESS", "KRS", "SFG", "ZEN",          # south fork salmon
    "IR1", "IR2", "IR3",                 # imnaha
    "LC1", "LC2",                        # lolo creek
    "MAR",                               # upper middle fork salmon
    "MR1", "UGR", "UGS", "WR1", "WR2"    # grande ronde / wallowa
  )) |>
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326
  )

# retain populations containing at least one juvenile release site
spsm_pop_trim = spsm_pop |>
  st_filter(
    st_transform(juv_mrr_sites, st_crs(spsm_pop)),
    .predicate = st_intersects
  ) |>
  select(TRT_POPID, POP_NAME, MPG) |>
  st_transform(4326)

# ensure all points are also WGS84
juv_mrr_sites = st_transform(juv_mrr_sites, 4326)
iptds_sf      = st_transform(iptds_sf, 4326)

# Lower Granite Dam
lgr_sf = tibble(
  site_code = "LGR",
  longitude = -117.4283,
  latitude = 46.66050
) |>
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326
  )

# ---------------------
# COLORS AND MAP EXTENT

# population colors by major population group (MPG)
chnk_mpg_col = colorFactor(
  palette = "Set2",
  domain = spsm_pop_trim$MPG
)

# site colors
mark_col  = "orangered"
iptds_col = "black"
lgr_col   = "gold"

# map extent includes populations, release sites, and IPTDS
map_bounds = st_bbox(c(
  st_geometry(spsm_pop_trim),
  st_geometry(juv_mrr_sites),
  st_geometry(iptds_sf)
))

# -------------
# BUILD LEAFLET
sr_juv_leaflet = leaflet() |>
  # base map
  addProviderTiles(providers$Esri.WorldTopoMap, group = "Topographic") |>
  # streams layer
  addPolylines(data = streams, color = "blue", weight = 1) |>
  # spring/summer Chinook population
  addPolygons(
    data = spsm_pop_trim,
    group = "Sp/Sum Chinook Pops",
    fillColor = ~chnk_mpg_col(MPG),
    fillOpacity = 0.3,
    stroke = TRUE,
    color = "black",
    weight = 1.5,
    opacity = 0.85,
    smoothFactor = 1,
  label = ~paste0(TRT_POPID, ": ", POP_NAME),
  popup = ~paste0(
    "<b>Spring/Summer Chinook Salmon</b><br>",
    "<b>Population:</b> ", POP_NAME, "<br>",
    "<b>Pop ID:</b> ", TRT_POPID, "<br>",
    "<b>MPG:</b> ", MPG
  ),
  highlightOptions = highlightOptions(
    color = "black",
    weight = 3,
    fillOpacity = 0.40,
    bringToFront = FALSE
  )) |>
  # iptds detection sites
  addCircleMarkers(
    data = iptds_sf,
    group = "IPTDS Sites",
    label = ~site_code,
    radius = 7,
    fillColor = iptds_col,
    fillOpacity = 0.95,
    color = "white",
    weight = 1.5,
    opacity = 1,
    stroke = TRUE,
    popup = ~paste0(
      "<b>IPTDS Detection Site</b><br>",
      "<b>Site:</b> ", site_code, "<br>",
      "<b>Name:</b> ", name, "<br>"
  )) |>
  # juvenile mark locations
  addCircleMarkers(
    data = juv_mrr_sites,
    group = "Juvenile Marking Sites",
    label = ~paste0(site_code, ": ", name),
    radius = 6,
    fillColor = mark_col,
    fillOpacity = 0.95,
    color = "white",
    weight = 1.5,
    opacity = 1,
    stroke = TRUE,
    popup = ~paste0(
      "<b>Juvenile Marking Site</b><br>",
      "<b>Site:</b> ", site_code, "<br>",
      "<b>Name:</b> ", name, "<br>",
      "<b>River km:</b> ", rkm, "<br>",
      "<b>Fish Tagged (all seasons):</b> ",
      scales::comma(n_tagged)
  )) |>
  # Lower Granite Dam
  addCircleMarkers(
    data = lgr_sf,
    group = "Lower Granite Dam",
    label = ~site_code,
    radius = 9,
    fillColor = lgr_col,  
    fillOpacity = 1,
    color = "black",
    weight = 2,
    opacity = 1,
    popup = "<b>Lower Granite Dam (LGR)</b>"
  ) |>
  # population MPG legend
  addLegend(
    position = "bottomright",
    pal = chnk_mpg_col,
    values = spsm_pop_trim$MPG,
    title = "MPG",
    group = "Sp/Sum Chinook Pops",
    opacity = 0.8
  ) |>
  # site legend
  addLegend(
    position = "bottomleft",
    colors = c(mark_col, iptds_col, lgr_col),
    labels = c(
      "Juvenile Marking Sites",
      "IPTDS Sites",
      "Lower Granite Dam)"
    ),
    title = "Monitoring Sites",
    opacity = 1
  ) |>
  # map title
  addControl(
    html = paste0(
      "<div style='background:rgba(255,255,255,0.95);",
      "padding:10px 14px;border-radius:5px;",
      "box-shadow:0 1px 5px #aaa;'>",
      "<b>2025 Juvenile Marking</b><br>",
      "<span style='font-size:12px;'>",
      "Snake River Spring/Summer Chinook",
      "</span></div>"
    ),
    position = "topleft"
  ) |>
  # layer controls
  addLayersControl(
    baseGroups = c("Topographic"),
    overlayGroups = c(
      "Sp/Sum Chinook Pops",
      "Juvenile Marking Sites",
      "IPTDS Sites"
    ),
    options = layersControlOptions(collapsed = FALSE)
  ) |>
  # set initial map extent
  fitBounds(
    lng1 = map_bounds[["xmin"]],
    lat1 = map_bounds[["ymin"]],
    lng2 = map_bounds[["xmax"]],
    lat2 = map_bounds[["ymax"]]
  )

# display map
sr_juv_leaflet

# ------------
# SAVE LEAFLET
saveWidget(
  sr_juv_leaflet,
  file = here("shiny", "leaflet", "sr_juv_marking_2025.html")
)

### END SCRIPT
