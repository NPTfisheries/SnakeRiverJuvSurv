# CREATE JUVENILE MARK GROUPS
# Inventory broadly enabled groups, evaluate YAML-selected sites, and retain
# one eligible release record per PIT tag in each qualifying group.

create_mark_groups = function(mark_events, config, lookups) {
  
  group_vars = c("species_run", "rear", "release_year", "method_group", "release_site", "season_life_stage")
  
  # classify events using the YAML species, method, and season definitions.
  raw_srr = dplyr::coalesce(as.character(mark_events$species_run_rear_type), "")
  species = rep(NA_character_, nrow(mark_events))
  for (i in seq_len(nrow(lookups$species_lookup))) {
    species[startsWith(raw_srr, lookups$species_lookup$prefix[i])] =
      lookups$species_lookup$species_run[i]
  }
  
  available_marks = mark_events |>
    dplyr::mutate(
      species_run = species,
      rear        = substr(raw_srr, 3, 3),
      release_md  = format(as.Date(release_date), "%m-%d")
    ) |>
    dplyr::filter(
      !is.na(species_run),
      rear %in% lookups$enabled_rear,
      is.null(lookups$release_years) | release_year %in% unlist(lookups$release_years, use.names = FALSE)
    ) |>
    dplyr::inner_join(
      lookups$method_lookup,
      by = c("capture_method", "rear"), 
      relationship = "many-to-one"
    ) |>
    dplyr::left_join(
      lookups$season_lookup,
      by = c("species_run", "release_md"), 
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(season_life_stage = tidyr::replace_na(season_life_stage, "unassigned"))
  
  # inventory all observed release sites BEFORE applying method-specific release sites.
  # n_start_tags counts valid PIT tags, not all release-event records.
  available_groups = available_marks |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      n_start_tags = dplyr::n_distinct(pit_tag[!is.na(pit_tag) & !stringr::str_trim(pit_tag) %in% c("", "..........")]),
      .groups = "drop"
    )
  
  allowed_sites = purrr::imap_dfr(
    config$method_group[lookups$enabled_methods],
    function(x, method) {
      if (identical(x$sites, "all")) {
        available_groups |>
          dplyr::filter(method_group == .env$method) |>
          dplyr::distinct(method_group, release_site)
      } else {
        tibble::tibble(method_group = method,
                       release_site = as.character(unlist(x$sites, use.names = FALSE)))
      }
    }
  ) |> dplyr::distinct()
  
  # flag events excluded from configured groups. exclusion counts may overlap.
  candidates = available_marks |>
    dplyr::semi_join(allowed_sites, by = c("method_group", "release_site")) |>
    dplyr::mutate(
      qc_invalid_pit = is.na(pit_tag) | stringr::str_trim(dplyr::coalesce(pit_tag, "")) %in% c("", ".........."),
      qc_missing_other = is.na(event_date) | is.na(release_date) | is.na(release_year) | is.na(release_site) |
        stringr::str_trim(dplyr::coalesce(release_site, "")) == "",
      qc_unassigned_season_life_stage = season_life_stage == "unassigned"
    )
  
  exclusion_cols = character()
  for (name in names(config$exclusions)) {
    rule    = config$exclusions[[name]]
    pattern = paste0("\\b(?:", paste(toupper(unlist(rule$codes)), collapse = "|"), ")\\b")
    applies = identical(rule$species, "all") | candidates$species_run %in% unlist(rule$species, use.names = FALSE)
    column  = paste0("qc_", name)
    candidates[[column]] = isTRUE(rule$enabled) & applies & stringr::str_detect(toupper(dplyr::coalesce(as.character(candidates$conditional_comments), "")), pattern)
    exclusion_cols = c(exclusion_cols, column)
  }
  
  qc_cols = c("qc_invalid_pit", "qc_missing_other", "qc_unassigned_season_life_stage", exclusion_cols)
  candidates = candidates |>
    dplyr::mutate(qc_exclude = dplyr::if_any(dplyr::all_of(qc_cols), identity))
  
  # summarize only selected sites. For non-selected sites, QC stays NA.
  selected_summary = candidates |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      n_invalid_pit           = sum(qc_invalid_pit),
      n_missing_or_unassigned = sum(qc_missing_other | qc_unassigned_season_life_stage),
      dplyr::across(dplyr::all_of(exclusion_cols), ~ sum(.x), .names = "n_{.col}"),
      n_excluded_events       = sum(qc_exclude),
      n_eligible_events       = sum(!qc_exclude),
      n_final_tags            = dplyr::n_distinct(pit_tag[!qc_exclude]),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      n_repeats_removed = n_eligible_events - n_final_tags,
      keep_group        = n_final_tags >= config$min_tags
    ) |>
    dplyr::select(-n_eligible_events) |>
    dplyr::rename_with(~ sub("^n_qc_", "n_", .x), dplyr::starts_with("n_qc_"))
  
  # a single preview table: all broadly available groups, with QC for selected sites.
  exclusion_counts = paste0("n_", names(config$exclusions))
  mark_group_summary = available_groups |>
    dplyr::left_join(selected_summary, by = group_vars, relationship = "one-to-one") |>
    dplyr::mutate(
      selected       = !is.na(keep_group),
      n_tags_removed = n_start_tags - n_final_tags
    ) |>
    dplyr::select(
      dplyr::all_of(group_vars), selected, n_start_tags, n_final_tags,
      n_tags_removed, keep_group, n_excluded_events, n_repeats_removed,
      n_invalid_pit, n_missing_or_unassigned, dplyr::all_of(exclusion_counts)
    ) |>
    dplyr::arrange(dplyr::across(dplyr::all_of(group_vars)))
  
  # repeat rule determines which eligible event represents each PIT in a group.
  primary    = if (endsWith(config$repeat_rule, "_event")) "event_date" else "release_date"
  secondary  = if (primary == "event_date") "release_date" else "event_date"
  descending = startsWith(config$repeat_rule, "last_")
  
  final_marks = candidates |>
    dplyr::filter(!qc_exclude) |>
    dplyr::arrange(
      if (descending) dplyr::desc(.data[[primary]]) else .data[[primary]],
      if (descending) dplyr::desc(.data[[secondary]]) else .data[[secondary]],
      file_name, sequence_number
    ) |>
    dplyr::distinct(dplyr::across(dplyr::all_of(c(group_vars, "pit_tag"))), .keep_all = TRUE) |>
    dplyr::semi_join(dplyr::filter(selected_summary, keep_group), by = group_vars) |>
    dplyr::select(-dplyr::all_of(c("release_md", "qc_exclude", qc_cols)))
  
  # each retained row corresponds to one eligible PIT tag within one group.
  stopifnot(nrow(final_marks) == sum(selected_summary$n_final_tags[selected_summary$keep_group]))
  
  list(mark_group_summary = mark_group_summary, final_marks = final_marks)
  
}

# END create_mark_groups()
