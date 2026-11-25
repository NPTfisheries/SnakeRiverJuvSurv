#' @title Compile Cached PTAGIS MRR Mark and Recapture Events
#'
#' @description
#' Compile previously cached PTAGIS MRR Mark and Recapture events into a single dataset for downstream mark-group creation.
#'
#' @details
#' This function discovers MRR caches created by \code{\link{cache_mrr_data}} and optionally filters them by project
#' code and PTAGIS file year. Setting \code{project_codes} or \code{cache_years} to NULL compiles all available caches.
#'
#' For each selected cache, session attributes are joined to event records using \code{file_name}. Only records with
#' \code{event_type} equal to "Mark" or "Recapture" are retained. Both event types represent potential release occasions.
#'
#' The original PTAGIS fields, including attached SPDV and PDV values, are retained. The function also adds \code{cache_project},
#' \code{cache_year}, \code{event_year}, and \code{release_year}. The latter two are derived from PTAGIS event and release dates.
#'
#' A source manifest is created which records the contributing caches, their MD5 checksums, and record counts. 
#' Original SPDV/PDV definitions remain in the source caches.
#'
#' @param cache_dir Parent directory containing cached PTAGIS MRR project/year RDS files.
#' @param project_codes Character vector of PTAGIS project codes to include. NULL includes all available projects.
#' @param cache_years Numeric vector of PTAGIS MRR file years to include. NULL includes all available years.
#'
#' @return A list containing:
#' \describe{
#'   \item{mark_events}{Tibble containing all selected Mark and Recapture events, joined with their session attributes
#'     and associated cache, event-year, and release-year fields.}
#'   \item{source_manifest}{Tibble containing one row per selected cache, with project code, file year, source path, MD5
#'     checksum, and counts of sessions, total events, and Mark/Recapture events.}
#' }
#'
#' @examples
#' \dontrun{
#' prepared <- compile_mrr_mark_events(
#'   cache_dir = "data/cache/ptagis/mrr",
#'   project_codes = c("CDR", "GAA"),
#'   cache_years = 2021:2025
#' )
#'
#' mark_events <- prepared$mark_events
#' source_manifest <- prepared$source_manifest
#' }
#'
#' @seealso \code{\link{cache_mrr_data}},
#'   \code{\link{create_mark_groups}}
#'
#' @author Mike Ackerman
#'
#' @export
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
