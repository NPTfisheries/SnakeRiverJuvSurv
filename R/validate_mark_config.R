# VALIDATE MARK-GROUP CONFIGURATION
# Check that YAML definitions are unambiguous; return lookup tables

validate_mark_config = function(config) {

  # check whether prep_id is valid
  if (!is.character(config$prep_id) || length(config$prep_id) != 1L || is.na(config$prep_id) || !nzchar(config$prep_id)) {
    stop("Set prep_id in the YAML to a prepared-data folder.")
  }
  
  # check whether repeat_rule is valid
  if (!is.character(config$repeat_rule) || length(config$repeat_rule) != 1L || is.na(config$repeat_rule) ||
      !config$repeat_rule %in% c("first_event", "last_event", "first_release", "last_release")) {
    stop("Invalid repeat_rule in YAML. Must be one of first_event, last_event, first_release, or last_release")
  }
  
  # check whether min_tag is valid
  if (length(config$min_tags) != 1L || !is.numeric(config$min_tags) ||
      is.na(config$min_tags) || !is.finite(config$min_tags) ||
      config$min_tags < 1 ||
      config$min_tags != trunc(config$min_tags)) {
    stop("min_tags must be a positive whole number.")
  }

  # check whether at least one species/run, method, and rear type is enabled
  enabled_species = names(config$species_run)[purrr::map_lgl(config$species_run, ~ isTRUE(.x$enabled))] 
  enabled_methods = names(config$method_group)[purrr::map_lgl(config$method_group, ~ isTRUE(.x$enabled))]
  enabled_rear    = names(config$rear)[purrr::map_lgl(config$rear, isTRUE)]
  
  if (!length(enabled_species) || !length(enabled_methods) || !length(enabled_rear)) {
    stop("Enable at least one species/run, method group, and rear type.")
  }

  # SRR prefixes must be unique and may not be nested across species.
  species_lookup = purrr::imap_dfr(
    config$species_run[enabled_species],
    ~ tibble::tibble(species_run = .y,
                     prefix = as.character(unlist(.x$srr_prefixes, use.names = FALSE)))
  )
  prefixes = species_lookup$prefix
  if (!length(prefixes) || anyNA(prefixes) || any(!grepl("^[A-Za-z0-9]+$", prefixes)) ||
      anyDuplicated(prefixes)) stop("Invalid or duplicate SRR prefixes in YAML.")
  for (i in seq_along(prefixes)) {
    if (any(startsWith(prefixes[-i], prefixes[i]))) {
      stop("Overlapping SRR prefixes in YAML.")
    }
  }

  # each capture_method + rear combination must belong to at most one group.
  method_lookup = purrr::imap_dfr(
    config$method_group[enabled_methods],
    function(x, method) {
      
      codes = as.character(unlist(x$codes, use.names = FALSE))
      rears = as.character(unlist(x$rear, use.names = FALSE))
      
      if (!length(codes) || !length(rears) || anyNA(codes) || anyNA(rears)) {
        stop("Define capture-method codes and rear types for: ", method)
      }
      if (!identical(x$sites, "all") && !length(x$sites)) {
        stop("Define sites or use 'all' for: ", method)
      }
      tidyr::crossing(capture_method = codes, rear = rears) |>
        dplyr::mutate(method_group = method)
    }
  )
  if (anyDuplicated(method_lookup[c("capture_method", "rear")])) {
    stop("Capture-method/rear definitions overlap in YAML.")
  }

  # a leap-year calendar supports whole months or recurring inclusive MM-DD ranges.
  calendar = tibble::tibble(day = seq(as.Date("2000-01-01"),
                                      as.Date("2000-12-31"), 
                                      by = "day")) |>
    dplyr::mutate(md = format(day, "%m-%d"), month = as.integer(format(day, "%m")))

  if (!all(enabled_species %in% names(config$season_life_stage))) {
    stop("Missing season_life_stage definitions for an enabled species.")
  }
  
  season_lookup = purrr::imap_dfr(
    config$season_life_stage[enabled_species],
    function(definitions, species) {
      purrr::imap_dfr(definitions, function(def, label) {
        if (identical(def, "all")) {
          include = rep(TRUE, nrow(calendar))
        } else if (is.list(def) && all(c("start", "end") %in% names(def))) {
          start = as.character(def$start)
          end = as.character(def$end)
          dates = suppressWarnings(as.Date(paste0("2000-", c(start, end))))
          if (length(start) != 1L || length(end) != 1L || anyNA(dates) ||
              !identical(format(dates, "%m-%d"), c(start, end))) {
            stop("Invalid MM-DD range for ", species, " / ", label)
          }
          include = if (start <= end) {
            calendar$md >= start & calendar$md <= end
          } else {
            calendar$md >= start | calendar$md <= end
          }
        } else {
          months = suppressWarnings(as.integer(unlist(def, use.names = FALSE)))
          if (!length(months) || anyNA(months) ||
              any(!months %in% 1:12) || anyDuplicated(months)) {
            stop("Invalid months for ", species, " / ", label)
          }
          include = calendar$month %in% months
        }
        tibble::tibble(species_run = species,
                       release_md = calendar$md[include],
                       season_life_stage = label)
      })
    }
  )
  if (anyDuplicated(season_lookup[c("species_run", "release_md")])) {
    stop("Overlapping season_life_stage definitions in YAML.")
  }

  # exclusion codes are matched as whole PTAGIS codes during group creation.
  for (name in names(config$exclusions)) {
    rule = config$exclusions[[name]]
    codes = unlist(rule$codes, use.names = FALSE)
    if (!length(codes) || anyNA(codes) ||
        any(!grepl("^[A-Za-z0-9]+$", as.character(codes)))) {
      stop("Invalid conditional-comment codes for: ", name)
    }
    if (!identical(rule$species, "all") &&
        !all(unlist(rule$species, use.names = FALSE) %in% names(config$species_run))) {
      stop("Unknown species in exclusion: ", name)
    }
  }

  list(
    species_lookup  = species_lookup,
    method_lookup   = method_lookup,
    season_lookup   = season_lookup,
    enabled_methods = enabled_methods,
    enabled_rear    = enabled_rear,
    release_years   = config$release_years
  )
}
