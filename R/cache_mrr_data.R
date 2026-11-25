#' @title Retrieve and Cache PTAGIS MRR Data
#'
#' @description
#' Retrieve and cache PTAGIS mark-recapture-recovery (MRR) data for a specified project code and file year.
#'
#' @details
#' The function uses \code{PITmodelR::get_mrr_files()} to identify submitted MRR files 
#' and \code{PITmodelR::get_batch_file_data()} to retrieve their contents.
#'
#' Data are cached separately by project and file year at \code{<cache_dir>/<project_code>/mrr_<project_code>_<year>.rds}.
#' Existing caches are reused unless \code{refresh = TRUE}, in which case the cache is replaced following successful retrieval.
#'
#' Each cached RDS contains the original \code{sessions}, \code{events}, and \code{pdv_map} tables, along with \code{file_inventory} 
#' and \code{retrieval} metadata. SPDV and PDV variables are attached but not mapped to descriptive column names.
#'
#' Failed retrievals produce a warning and return a status indicating failure. Previously cached data are not intentionally 
#' removed when retrieval fails.
#'
#' @param project_code Three-character PTAGIS MRR project code.
#' @param year Four-digit PTAGIS MRR file year.
#' @param cache_dir Parent directory for cached MRR data.
#' @param refresh Logical. If TRUE, retrieve and replace an existing cache. Defaults to FALSE.
#' @param page_size Number of file records requested per API page. Defaults to 1000.
#'
#' @return A one-row tibble containing:
#' \describe{
#'   \item{project_code}{PTAGIS project code.}
#'   \item{year}{PTAGIS MRR file year.}
#'   \item{status}{One of "cached", "downloaded", "refreshed", "no_files", or "failed".}
#'   \item{n_files}{Number of MRR files retrieved.}
#'   \item{n_sessions}{Number of session records retrieved.}
#'   \item{n_events}{Number of event records retrieved.}
#'   \item{cache_path}{Path to the cached RDS file.}
#'   \item{attempted_dttm}{Formatted retrieval-attempt date and time.}
#'   \item{retrieved_dttm}{Formatted successful-retrieval date and time.}
#'   \item{error}{Error message if retrieval failed.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- cache_mrr_data(
#'   project_code = "CDR",
#'   year = 2025,
#'   cache_dir = "data/cache/ptagis/mrr"
#' )
#'
#' mrr <- readRDS(result$cache_path)
#' }
#'
#' @seealso \code{PITmodelR::get_mrr_files()},
#'   \code{PITmodelR::get_batch_file_data()},
#'   \code{\link{compile_mrr_mark_events}}
#'
#' @author Mike Ackerman
#'
#' @export
cache_mrr_data <- function(project_code,
                           year,
                           cache_dir,
                           refresh = FALSE,
                           page_size = 1000L) {
  
  # validate inputs
  if (!is.character(project_code) || length(project_code) != 1L || is.na(project_code) || nchar(trimws(project_code)) != 3L) {
    stop("`project_code` must be a three-character string.")
  }
  
  if (!is.numeric(year) || length(year) != 1L || is.na(year) || !is.finite(year) || year != trunc(year) || year < 1000 || year > 9999) {
    stop("`year` must be a four-digit integer.")
  }
  
  if (!is.character(cache_dir) || length(cache_dir) != 1L || is.na(cache_dir) || !nzchar(trimws(cache_dir))) {
    stop("`cache_dir` must be a nonempty directory path.")
  }
  
  if (!is.logical(refresh) || length(refresh) != 1L || is.na(refresh)) {
    stop("`refresh` must be TRUE or FALSE.")
  }
  
  project_code <- toupper(trimws(project_code))
  year         <- as.integer(year)
  
  # cache location
  cache_path <- file.path(
    cache_dir,
    project_code,
    paste0("mrr_", project_code, "_", year, ".rds")
  )
  
  had_cache <- file.exists(cache_path)
  
  attempted_dttm <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z",, tz = "America/Denver")
  
  # result summary helper
  make_result <- function(status,
                          n_files           = NA_integer_,
                          n_sessions        = NA_integer_,
                          n_events          = NA_integer_,
                          retrieved_dttm    = NA_character_,
                          error             = NA_character_) {
    
    tibble::tibble(
      project_code = project_code,
      year = year,
      status = status,
      n_files = n_files,
      n_sessions = n_sessions,
      n_events = n_events,
      attempted_dttm = attempted_dttm,
      retrieved_dttm = retrieved_dttm,
      cache_path     = cache_path,
      error          = error
    )
  }
  
  # skip existing cache unless refresh = TRUE
  if (had_cache && !refresh) {
    message("Cached MRR data already exists: ", project_code, " / ", year)
    return(make_result("cached"))
  }
  
  # retrieve and cache
  tryCatch({
    
    message("Retrieving MRR data: ", project_code, " / ", year)
    
    # retrieve file inventory
    files <- PITmodelR::get_mrr_files(
      code      = project_code,
      year      = year,
      page_size = page_size,
      all_pages = TRUE
    )
    
    if (!is.data.frame(files)) {
      stop("Unexpected MRR file inventory.")
    }
    
    if (nrow(files) == 0L) {
      return(make_result("no_files", n_files = 0L))
    }
    
    if (!"name" %in% names(files) || anyNA(files$name) || any(!nzchar(trimws(files$name))) || anyDuplicated(files$name)) {
      stop("Missing or duplicate filenames in MRR inventory.")
    }
    
    # retrieve MRR records
    mrr <- PITmodelR::get_batch_file_data(
      filenames = files$name,
      pdvs      = "attach",
      map_pdvs  = FALSE
    )
    
    # basic retrieval QC
    if (!is.list(mrr) || !is.data.frame(mrr$sessions) || !is.data.frame(mrr$events)) {
      stop("Unexpected MRR data structure.")
    }
    
    if (!"file_name" %in% names(mrr$sessions) || !"file_name" %in% names(mrr$events)) {
      stop("Missing file_name in MRR records.")
    }
    
    if (nrow(mrr$sessions) != nrow(files) || anyDuplicated(mrr$sessions$file_name) || !setequal(mrr$sessions$file_name, files$name)) {
      stop("Session records do not match file inventory.")
    }
    
    if (!all(mrr$events$file_name %in% files$name)) {
      stop("Events contain unexpected filenames.")
    }
    
    # preserve PDV map, including when no PDVs are returned
    if (is.null(mrr$pdv_map)) {
      mrr$pdv_map <- tibble::tibble()
    }
    
    if (!is.data.frame(mrr$pdv_map)) {
      stop("Unexpected PDV map structure.")
    }
    
    if (nrow(mrr$pdv_map) > 0L && (!"file_name" %in% names(mrr$pdv_map) || !all(mrr$pdv_map$file_name %in% files$name))) {
      stop("PDV map contains unexpected filenames.")
    }
    
    # add retrieval metadata
    retrieved_dttm <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z",, tz = "America/Denver")
    
    mrr$file_inventory <- files
    
    mrr$retrieval <- list(
      project_code      = project_code,
      year              = year,
      retrieved_dttm    = retrieved_dttm,
      PITmodelR_version = as.character(utils::packageVersion("PITmodelR")),
      pdvs              = "attach",
      map_pdvs          = FALSE
    )
    
    # save cache
    dir.create(
      dirname(cache_path),
      recursive    = TRUE,
      showWarnings = FALSE
    )
    
    saveRDS(mrr, cache_path)
    
    message("Saved: ", cache_path)
    
    # return retrieval summary
    make_result(
      status            = if (had_cache) "refreshed" else "downloaded",
      n_files           = nrow(files),
      n_sessions        = nrow(mrr$sessions),
      n_events          = nrow(mrr$events),
      retrieved_dttm    = retrieved_dttm
    )
    
  }, error = function(e) {
    
    warning("MRR retrieval failed for ", project_code, " / ", year, ": ", conditionMessage(e), call. = FALSE)
    
    make_result(
      status = "failed",
      error = conditionMessage(e)
    )
  })
} # END cache_mrr_data()
