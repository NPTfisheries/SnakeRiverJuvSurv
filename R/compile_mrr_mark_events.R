# COMPILE CACHED PTAGIS MRR MARK EVENTS
#
# Read selected project/year caches, join sessions to events within each cache,
# and retain all Mark and Recapture records for downstream mark-group creation.

compile_mrr_mark_events = function(cache_dir, 
                                   project_codes = NULL, 
                                   cache_years = NULL) {

  # discover all available caches; NULL retains all; project_codes and cache_years filtered below
  cache_files = list.files(
    cache_dir,
    pattern    = "^mrr_[A-Z0-9]{3}_[0-9]{4}\\.rds$",
    recursive  = TRUE,
    full.names = TRUE
  ) |> sort()

  if (length(cache_files) == 0L) stop("No MRR caches found in: ", cache_dir)

  cache_index = tibble::tibble(
    cache_project = stringr::str_match(basename(cache_files), "^mrr_([A-Z0-9]{3})_([0-9]{4})\\.rds$")[, 2],
    cache_year    = as.integer(stringr::str_match(basename(cache_files), "^mrr_([A-Z0-9]{3})_([0-9]{4})\\.rds$")[, 3]),
    cache_file    = normalizePath(cache_files)
  ) |>
    # filter to project_codes and cache_years if not NULL
    dplyr::filter(
      is.null(project_codes) | cache_project %in% project_codes,
      is.null(cache_years)   | cache_year %in% cache_years
    )

  if (nrow(cache_index) == 0L) stop("No MRR caches match the selected projects and years.")

  # warn only when an explicitly requested project/year is not found
  if (!is.null(project_codes)) {
    missing = setdiff(project_codes, cache_index$cache_project)
    if (length(missing)) warning("No MRR caches found for project(s): ", paste(missing, collapse = ", "))
  }
  if (!is.null(cache_years)) {
    missing = setdiff(cache_years, cache_index$cache_year)
    if (length(missing)) warning("No MRR caches found for year(s): ", paste(missing, collapse = ", "))
  }
  
  # md5 can be used if we ever need to track whether the contents of a cached file has changed
  cache_index = dplyr::mutate(cache_index, md5 = unname(tools::md5sum(cache_file)))

  cache_results = purrr::map(seq_len(nrow(cache_index)), function(i) {
    
    mrr  = readRDS(cache_index$cache_file[i])

    sessions = tibble::as_tibble(mrr$sessions)
    events   = tibble::as_tibble(mrr$events)
    
    mark_events = events |>
      dplyr::filter(event_type %in% c("Mark", "Recapture")) |>
      dplyr::left_join(
        sessions, 
        by           = "file_name", 
        relationship = "many-to-one",
        suffix       = c("", "_session")
      ) |>
      dplyr::mutate(
        cache_project = cache_index$cache_project[i],
        cache_year    = cache_index$cache_year[i],
        event_year    = lubridate::year(event_date),
        release_year  = lubridate::year(release_date)
      )

    list(
      mark_events         = mark_events,
      n_sessions          = nrow(sessions),
      n_events            = nrow(events),
      n_mark_recap_events = nrow(mark_events)
    )
  })

  mark_events = purrr::list_rbind(purrr::map(cache_results, "mark_events"))
  if (nrow(mark_events) == 0L) stop("No Mark or Recapture events found.")

  source_manifest = cache_index |>
    dplyr::mutate(
      n_sessions          = purrr::map_int(cache_results, "n_sessions"),
      n_events            = purrr::map_int(cache_results, "n_events"),
      n_mark_recap_events = purrr::map_int(cache_results, "n_mark_recap_events")
    )

  stopifnot(nrow(mark_events) == sum(source_manifest$n_mark_recap_events))

  list(
    mark_events     = mark_events, 
    source_manifest = source_manifest
  )
}

# END compile_mrr_mark_events()
