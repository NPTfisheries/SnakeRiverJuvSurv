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

check_mark_groups_tag_yr <- function(activity = c('parr', 'rst', 'hatchery'), event_type = 'Mark', tag_year = NULL) {
  
  activity = match.arg(activity)
  
  stopifnot(event_type == 'Mark')#{'event_type must be equal to "Mark"'}
  
  stopifnot(!is.null(tag_year)|!is.integer(tag_year))#{'tag_year must be a postive integer'}
  
  file_path <<- paste0('./output/events_', activity, "/tag_yr", tag_yr, "/", event_type)
  
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
          p_formula = ~stratum:time,
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
      cjs_p$interval <- cjs_p$interval + 1 #TEMPORARY FIX
      
      mscjs_p <- fit_mscjs$p 
      mscjs_p$interval <- mscjs_p$interval + 1 #TEMPORARY FIX
 
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


run_cjs_model_set_multiple <- function(ch_long_list,
                                       min_detected_occasions = 3) {
  
  safe_build <- safely(build_capture_histories)
  safe_cjs <- safely(fit_marked_cjs_multiple)
  safe_mscjs <- safely(fit_marked_mscjs_multiple)
  
  count_detected_occasions <- function(ch_data) {
    if (is.null(ch_data) || !"ch" %in% names(ch_data) || !nrow(ch_data)) return(0L)
    
    ch_mat <- do.call(rbind, strsplit(as.character(ch_data$ch), ""))
    sum(colSums(ch_mat == "1" | ch_mat == "2") > 0)
  }
  
  extract_error <- function(x) {
    if (is.null(x$error)) NA_character_ else x$error$message
  }
  
  ch_list <- ch_long_list %>%
    mutate(
      ch = map2(
        data,
        obs_loc,
        ~ build_capture_histories(
          tag_history = .x,
          locs_def = .y,
          site_col = "obs_site",
          tag_col = "tag_code",
          time_col = "first_obs",
          censor_col = "last_censored",
          covariate_cols = NULL
        )
      )
    )
  
  marray_list <- ch_list %>%
    mutate(m_array = map(
      ch, ~(marray = .x$m_array))) %>%
    select(srr, release_site, m_array) %>%
    rename(release_group = release_site)
  
  site_mapping <- ch_list %>%
    mutate(site_map = map(ch, "mapping")) %>%
    select(release_site, site_map) %>%
    rename(release_group = release_site) %>%
    unnest(site_map) %>%
    ungroup() %>%
    select(release_group, occasion, occ_idx) %>%
    distinct() 
  
  max_occ <- max(site_mapping$occ_idx)
  
  #adjust occasion indexing for release groups with a lower number of occasions
  # this is important for sharing information for sites downstream (i.e., LGR and Down)
  site_mapping <- site_mapping %>%
    group_by(release_group) %>%
    mutate(occ_idx = occ_idx + (max_occ - max(occ_idx))) %>%
    ungroup() %>%
    rename(occ = occ_idx, 
           detect_site = occasion) 
  
  multiple_ch_list <- ch_list %>%
    group_by(srr, release_season) %>% 
    #combine release sites together, keep separate by srr and release season
    summarise(ch_data = list(map2_dfr(ch, release_site, ~ mutate(.x$ch_data, release_site = .y))), 
              .groups = "drop") %>%
    #determine length of longest ch by group
    mutate(ch_data = map(ch_data, ~ {max_len <- max(nchar(.x$ch)) 
    .x %>%
      #add leading 0s to shorter ch, all release groups end with detections at LGR and downstream
      mutate(ch = str_pad(ch, width = max_len, side = "left", pad = "0"))}), 
    ch_data = map(ch_data, ~ mutate(.x, release_group = as.factor(.x$release_site))), 
    ch_data = map(ch_data, as.data.frame)
    )
  
  model_tbl <- multiple_ch_list %>%
    mutate(
      error_build = map_chr(ch_data, extract_error),
      
      n_detected_occasions = map_int(
        ch_data,
        ~ count_detected_occasions(.x)
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
      
      cjs_data = map2(ch_data, run_models, ~ {
        if (!.y) return(NULL)
        x <- .x
        x$ch <- gsub("2", "1", x$ch)
        x
      }),
      
      ms_data = map2(ch_data, run_models, ~ {
        if (!.y) return(NULL)
        build_multistate_histories(.x)
      }),
      
      fit_cjs_safe = map2(
        cjs_data,
        run_models,
        ~ if (.y) safe_cjs(.x, phi_formula = ~time*release_group, p_formula = ~time*detect_site, 
                           site_mapping = site_mapping) else list(result = NULL, error = NULL)
      ),
      
      fit_mscjs_safe = map2(
        ms_data,
        run_models,
        ~ if (.y) safe_mscjs(
          .x,
          site_mapping = site_mapping,
          s_formula = ~time*release_group,
          p_formula = ~stratum*time*detect_site,
          psi_formula = ~ -1 + stratum:tostratum
        ) else list(result = NULL, error = NULL)
      ),
      
      fit_cjs = map(fit_cjs_safe, "result"),
      fit_mscjs = map(fit_mscjs_safe, "result"),
      # fit_marray = map(fit_marray_safe, "result"),
      
      error_cjs = map_chr(fit_cjs_safe, extract_error),
      error_mscjs = map_chr(fit_mscjs_safe, extract_error),
    )
  
  fit_marray <- fit_marray_cjs_multiple(marray_list)
  
  estimate_df <- model_tbl %>%
    select(
      srr, release_season, #obs_loc,
      n_detected_occasions,
      error_build, error_precheck, error_cjs,
      error_mscjs,
      #error_marray,
      fit_cjs, fit_mscjs
      #fit_marray
    ) %>%
    pmap_dfr(function(srr,
                      release_season,
                      #obs_loc,
                      n_detected_occasions,
                      error_build,
                      error_precheck,
                      error_cjs,
                      error_mscjs,
                      #error_marray,
                      fit_cjs,
                      fit_mscjs
                      #fit_marray
    ) {
      
      #loc_names <- names(obs_loc)
      #n_sites <- length(obs_loc)
      n_reaches <- n_detected_occasions - 1
      n_release_group <- length(unique(site_mapping$release_group))
      sites <- site_mapping %>%
        filter(occ > 1) %>%
        select(detect_site) %>%
        distinct()
      
      group_cols <- tibble(
        srr = srr,
        #release_site = release_site,
        release_season = release_season,
        n_detected_occasions = n_detected_occasions,
        error_build = error_build,
        error_precheck = error_precheck,
        error_cjs = error_cjs,
        error_mscjs = error_mscjs,
        #error_marray = error_marray
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
                   mscjs_ucl = NA_real_
                   # marray_est = NA_real_,
                   # marray_lcl = NA_real_,
                   # marray_ucl = NA_real_
                 ))
      }
      
      cjs_phi <- fit_cjs$phi
      mscjs_phi <- fit_mscjs$phi
      marray_phi <- fit_marray$phi %>%
        mutate(metric = "survival") %>%
        rename(marray_est = phi, 
               marray_lcl = lcl, 
               marray_ucl = ucl)
      
      cjs_cum <- fit_cjs$cum_phi
      mscjs_cum <- fit_mscjs$cum_phi
      #marray_cum <- fit_marray$cum_phi
      
      cjs_p <- fit_cjs$p
      #cjs_p$interval <- cjs_p$interval
      
      mscjs_p <- fit_mscjs$p 
      #mscjs_p$interval <- mscjs_p$interval + 1 #TEMPORARY FIX
      
      marray_p <- fit_marray$p %>%
        mutate(metric = "detection") %>%
        rename(marray_est = p, 
               marray_lcl = lcl, 
               marray_ucl = ucl, 
               to_site = site)
      
      surv_df <- cjs_phi %>%
        mutate(metric = "survival") %>%
        left_join(site_mapping %>% mutate(interval = as.factor(occ))) %>%
        mutate(occ = occ + 1) %>%
        rename(from_site = detect_site) %>%
        left_join(site_mapping, by = c("release_group", "occ")) %>%
        rename(to_site = detect_site) %>%
        select(-occ, -se) %>%
        rename(cjs_est = estimate,
               cjs_lcl = lcl,
               cjs_ucl = ucl) %>%
        left_join(mscjs_phi, by = c("interval", "release_group"), copy = T) %>%
        rename(mscjs_est = estimate,
               mscjs_lcl = lcl,
               mscjs_ucl = ucl) %>%
        select(-se) %>%
        left_join(marray_phi, by = c("release_group", "to_site", "from_site", "metric"))
      
      #left_join(make_est_df(marray_phi, "marray"), by = "interval")
      
      cum_df <- cjs_cum %>%
        mutate(metric = "cumulative",
               occ = as.numeric(interval) + 1) %>%
        left_join(site_mapping) %>%
        rename(to_site = detect_site) %>%
        mutate(from_site = release_group) %>%
        select(-occ) %>%
        rename(cjs_est = estimate,
               cjs_lcl = lcl,
               cjs_ucl = ucl) %>%
        left_join(mscjs_phi, by = c("interval", "release_group"), copy = T) %>%
        rename(mscjs_est = estimate,
               mscjs_lcl = lcl,
               mscjs_ucl = ucl)
      
      #left_join(make_est_df(marray_cum, "marray"), by = "interval")
      
      p_df <- cjs_p %>%
        mutate(metric = "detection") %>%
        rename(cjs_est = estimate,
               cjs_lcl = lcl,
               cjs_ucl = ucl) %>%
        select(-se) %>%
        left_join(mscjs_p, by = c("interval", "detect_site"), copy = T) %>%
        rename(mscjs_est = estimate,
               mscjs_lcl = lcl,
               mscjs_ucl = ucl,
               from_site = detect_site) %>%
        mutate(to_site = from_site) %>%
        left_join(marray_p, by = c("to_site", "metric"))
      
      #left_join(make_est_df(marray_p, "marray"), by = "interval")
      
      bind_rows(surv_df, cum_df, p_df) %>%
        bind_cols(group_cols[rep(1, nrow(.)), ]) #%>%
      #relocate(srr, release_site, release_season, n_detected_occasions)
    })
  
  list(
    model_tbl = model_tbl,
    estimates = estimate_df
  )
}