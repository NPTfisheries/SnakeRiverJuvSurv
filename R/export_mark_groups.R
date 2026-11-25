# EXPORT JUVENILE MARK GROUPS
# Save a snapshot of the configuration and group results, plus one or more
# headerless PIT-tag lists per qualifying group. Keep mark data in one RDS.

export_mark_groups = function(results, 
                              config, 
                              config_file,
                              output_dir = here::here("data", "processed", "mark_groups"),
                              max_tags_per_file = 15000L) {

  group_vars = c("species_run", "rear", "method_group", "season_life_stage", "release_site", "release_year")

  # only export qualifying mark groups; check that their records match the summary.
  marks  = results$final_marks
  groups = results$mark_group_summary |>
    dplyr::filter(keep_group %in% TRUE) |>
    dplyr::select(dplyr::all_of(group_vars), n_final_tags)

  if (!nrow(marks)) stop("No qualifying mark groups to export.")

  counts = marks |>
    dplyr::count(dplyr::across(dplyr::all_of(group_vars)), name = "n_export_tags")
  
  check = groups |>
    dplyr::left_join(
      counts, by = group_vars, 
      relationship = "one-to-one"
    )
  if (nrow(counts) != nrow(groups) || anyNA(check$n_export_tags) || any(check$n_export_tags != check$n_final_tags)) {
    stop("Final marks do not match the qualifying group summary.")
  }

  # create a mark group ID
  make_group_id = function(x) {
    id = paste(x$species_run, x$rear, x$method_group, x$season_life_stage, x$release_site, x$release_year, sep = "_")
    gsub("[^A-Za-z0-9_-]", "-", id)
  }
  groups$group_id = make_group_id(groups)
  if (anyDuplicated(groups$group_id)) stop("Group IDs collide after filename cleanup.")

  # create export directory
  export_id  = paste0("export_", format(Sys.time(), "%Y%m%d_%H%M"))
  export_dir = file.path(output_dir, config$analysis_id, export_id)
  dir.create(file.path(export_dir, "tag_lists"), recursive = TRUE)

  # write tag lists; each file has one tag per line, no header; split only for large groups
  tags_by_group = split(as.character(marks$pit_tag), make_group_id(marks))
  tag_files = purrr::map_dfr(seq_len(nrow(groups)), function(i) {
    group_id = groups$group_id[i]
    tags = sort(tags_by_group[[group_id]])
    n_parts = ceiling(length(tags) / max_tags_per_file)

    purrr::map_dfr(seq_len(n_parts), function(part) {
      start = (part - 1L) * max_tags_per_file + 1L
      end = min(part * max_tags_per_file, length(tags))
      part_tags = tags[start:end]
      filename = if (n_parts == 1L) {
        paste0(group_id, "_tag_list.txt")
      } else {
        sprintf("%s_part%02d_tag_list.txt", group_id, part)
      }
      writeLines(part_tags, file.path(export_dir, "tag_lists", filename))
      tibble::tibble(group_id = group_id, part = part,
                     file = paste0("tag_lists/", filename),
                     n_tags_file = length(part_tags))
    })
  })

  # save copy of YAML config, mark group summary, tag list manifest, and one consolidated mark dataset
  if (!file.copy(config_file, file.path(export_dir, "config.yml"))) {
    stop("Could not copy YAML configuration to export folder.")
  }
  readr::write_csv(results$mark_group_summary, file.path(export_dir, "mark_group_summary.csv"))
  saveRDS(marks, file.path(export_dir, "final_marks.rds"))
  readr::write_csv(
    groups |>
      dplyr::select(-n_final_tags) |>
      dplyr::left_join(tag_files, by = "group_id", relationship = "one-to-many") |>
      dplyr::relocate(group_id),
    file.path(export_dir, "tag_list_manifest.csv")
  )

  message("Exported ", nrow(groups), " mark groups (", nrow(tag_files), " tag files; ", nrow(marks), " group-level marks): ", export_dir)
  invisible(export_dir)
}

# END export_mark_groups()
