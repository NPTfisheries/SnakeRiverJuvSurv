# helper function
check_mark_groups <- function(activity = c('parr', 'rst', 'hatchery'), event_type = 'Mark', migratory_year = NULL) {
  
  activity = match.arg(activity)
  
  stopifnot(event_type == 'Mark')#{'event_type must be equal to "Mark"'}
  
  stopifnot(!is.null(migratory_year)|!is.integer(migratory_year))#{'migratory_year must be a postive integer'}
  
  file_path <<- paste0('./output/events_', activity, "/MY", mig_yr, "/", event_type)
  
  # find available mark groups
  mark_groups <- list.files(
    path = file_path,
    pattern = "Mark\\.rds$"
  )
  
  # expected matching INT files
  int_files <- sub("Mark\\.rds$", "INT.csv", mark_groups)
  
  cth_files <- sub("Mark\\.rds$", "CTH.csv", mark_groups)
  
  # check existence
  int_exists <- file.exists(
    file.path(file_path, int_files)
  )
  
  cth_exists <- file.exists(
    file.path(file_path, cth_files)
  )
  
  # return summary table
  data.frame(
    mark_file = mark_groups,
    int_file  = int_files,
    int_exists    = int_exists,
    cth_file = cth_files,
    cth_exists = cth_exists,
    stringsAsFactors = FALSE
  )
}


# helper runs all three models and compiles data
run_cjs_model_set <- function(ch_long_list,
                              min_detected_occasions = 3) {
  
  safe_build <- safely(build_capture_histories)
  safe_cjs <- safely(fit_marked_cjs)
  safe_mscjs <- safely(fit_marked_mscjs)
  safe_marray <- safely(fit_marray_cjs)
  
  make_est_df <- function(x, prefix) {
    if (is.null(x) || !is.data.frame(x) || !nrow(x)) {
      return(tibble(
        interval = integer(),
        "{prefix}_est" := numeric(),
        "{prefix}_lcl" := numeric(),
        "{prefix}_ucl" := numeric()
      ))
    }
    
    x %>%
      select(interval, estimate, lcl, ucl) %>%
      rename(
        "{prefix}_est" := estimate,
        "{prefix}_lcl" := lcl,
        "{prefix}_ucl" := ucl
      )
  }
  
  count_detected_occasions <- function(ch_data) {
    if (is.null(ch_data) || !"ch" %in% names(ch_data) || !nrow(ch_data)) return(0L)
    
    ch_mat <- do.call(rbind, strsplit(as.character(ch_data$ch), ""))
    sum(colSums(ch_mat == "1" | ch_mat == "2") > 0)
  }
  
  extract_error <- function(x) {
    if (is.null(x$error)) NA_character_ else x$error$message
  }
  
  model_tbl <- ch_long_list %>%
    mutate(
      ch_safe = map2(
        data,
        obs_loc,
        ~ safe_build(
          tag_history = .x,
          locs_def = .y,
          site_col = "obs_site",
          tag_col = "tag_code",
          time_col = "first_obs",
          censor_col = "last_censored",
          covariate_cols = NULL
        )
      ),
      ch_list = map(ch_safe, "result"),
      error_build = map_chr(ch_safe, extract_error),
      
      n_detected_occasions = map_int(
        ch_list,
        ~ count_detected_occasions(.x$ch_data)
      ),
      
      error_precheck = if_else(
        is.na(error_build) & n_detected_occasions < min_detected_occasions,
        paste0(
          "Only ", n_detected_occasions,
          " detected occasion(s); need at least ",
          min_detected_occasions, "."
        ),
        NA_character_
      ),
      
      run_models = is.na(error_build) & is.na(error_precheck),
      
      cjs_data = map2(ch_list, run_models, ~ {
        if (!.y) return(NULL)
        x <- .x$ch_data
        x$ch <- gsub("2", "1", x$ch)
        x
      }),
      
      ms_data = map2(ch_list, run_models, ~ {
        if (!.y) return(NULL)
        build_multistate_histories(.x$ch_data)
      }),
      
      m_array = map2(ch_list, run_models, ~ {
        if (!.y) return(NULL)
        .x$m_array
      }),
      
      fit_cjs_safe = map2(
        cjs_data,
        run_models,
        ~ if (.y) safe_cjs(.x, phi_formula = ~time, p_formula = ~time) else list(result = NULL, error = NULL)
      ),
      
      fit_mscjs_safe = map2(
        ms_data,
        run_models,
        ~ if (.y) safe_mscjs(
          .x,
          s_formula = ~time,
          p_formula = ~time,
          psi_formula = ~ -1 + stratum:tostratum
        ) else list(result = NULL, error = NULL)
      ),
      
      fit_marray_safe = map2(
        m_array,
        run_models,
        ~ if (.y) safe_marray(.x) else list(result = NULL, error = NULL)
      ),
      
      fit_cjs = map(fit_cjs_safe, "result"),
      fit_mscjs = map(fit_mscjs_safe, "result"),
      fit_marray = map(fit_marray_safe, "result"),
      
      error_cjs = map_chr(fit_cjs_safe, extract_error),
      error_mscjs = map_chr(fit_mscjs_safe, extract_error),
      error_marray = map_chr(fit_marray_safe, extract_error)
    )
  
  estimate_df <- model_tbl %>%
    select(
      srr, release_site, release_season, obs_loc,
      n_detected_occasions,
      error_build, error_precheck, error_cjs, error_mscjs, error_marray,
      fit_cjs, fit_mscjs, fit_marray
    ) %>%
    pmap_dfr(function(srr,
                      release_site,
                      release_season,
                      obs_loc,
                      n_detected_occasions,
                      error_build,
                      error_precheck,
                      error_cjs,
                      error_mscjs,
                      error_marray,
                      fit_cjs,
                      fit_mscjs,
                      fit_marray) {
      
      loc_names <- names(obs_loc)
      n_sites <- length(obs_loc)
      n_reaches <- n_sites - 1
      
      group_cols <- tibble(
        srr = srr,
        release_site = release_site,
        release_season = release_season,
        n_detected_occasions = n_detected_occasions,
        error_build = error_build,
        error_precheck = error_precheck,
        error_cjs = error_cjs,
        error_mscjs = error_mscjs,
        error_marray = error_marray
      )
      
      if (!is.na(error_build) || !is.na(error_precheck)) {
        return(group_cols %>%
                 mutate(
                   metric = "error",
                   interval = NA_integer_,
                   from_site = NA_character_,
                   to_site = NA_character_,
                   cjs_est = NA_real_,
                   cjs_lcl = NA_real_,
                   cjs_ucl = NA_real_,
                   mscjs_est = NA_real_,
                   mscjs_lcl = NA_real_,
                   mscjs_ucl = NA_real_,
                   marray_est = NA_real_,
                   marray_lcl = NA_real_,
                   marray_ucl = NA_real_
                 ))
      }
      
      cjs_phi <- fit_cjs$phi
      mscjs_phi <- fit_mscjs$phi
      marray_phi <- fit_marray$phi
      
      cjs_cum <- fit_cjs$cum_phi
      mscjs_cum <- fit_mscjs$cum_phi
      marray_cum <- fit_marray$cum_phi
      
      cjs_p <- fit_cjs$p
      mscjs_p <- fit_mscjs$p
      marray_p <- fit_marray$p
      
      surv_df <- tibble(
        metric = "survival",
        interval = seq_len(n_reaches),
        from_site = loc_names[seq_len(n_reaches)],
        to_site = loc_names[seq(2, n_sites)]
      ) %>%
        left_join(make_est_df(cjs_phi, "cjs"), by = "interval") %>%
        left_join(make_est_df(mscjs_phi, "mscjs"), by = "interval") %>%
        left_join(make_est_df(marray_phi, "marray"), by = "interval")
      
      cum_df <- tibble(
        metric = "cumulative",
        interval = seq_len(n_reaches),
        from_site = release_site,
        to_site = loc_names[seq(2, n_sites)]
      ) %>%
        left_join(make_est_df(cjs_cum, "cjs"), by = "interval") %>%
        left_join(make_est_df(mscjs_cum, "mscjs"), by = "interval") %>%
        left_join(make_est_df(marray_cum, "marray"), by = "interval")
      
      p_df <- tibble(
        metric = "detection",
        interval = seq_len(n_sites),
        from_site = loc_names,
        to_site = loc_names
      ) %>%
        left_join(make_est_df(cjs_p, "cjs"), by = "interval") %>%
        left_join(make_est_df(mscjs_p, "mscjs"), by = "interval") %>%
        left_join(make_est_df(marray_p, "marray"), by = "interval")
      
      bind_rows(surv_df, cum_df, p_df) %>%
        bind_cols(group_cols[rep(1, nrow(.)), ]) %>%
        relocate(srr, release_site, release_season, n_detected_occasions)
    })
  
  list(
    model_tbl = model_tbl,
    estimates = estimate_df
  )
}
