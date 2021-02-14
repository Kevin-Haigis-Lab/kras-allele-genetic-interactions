# A system for extracting data in plots and saving as "Source Data".

SOURCE_DATA_BASE_DIR <- here::here("paper", "source-data")

if (!dir.exists(SOURCE_DATA_BASE_DIR)) {
  dir.create(SOURCE_DATA_BASE_DIR)
}


#### ---- Data extraction ---- ####

pull_original_plot_data <- function(p) {
  as_tibble(p$data)
}


pull_wrapped_plot_data <- function(wp, annotation_col_name, annotations) {
  df <- map_dfr(seq(1, length(wp$patches$plots)), function(i) {
    a <- annotations[[i]]
    wp$patches$plots[[i]]$data %>%
      as_tibble() %>%
      mutate({{ annotation_col_name }} := !!a)
  })

  a <- annotations[[length(annotations)]]
  last_df <- pull_original_plot_data(wp) %>%
    as_tibble() %>%
    mutate({{ annotation_col_name }} := !!a)

  bind_rows(df, last_df)
}


remove_ggraph_columns <- function(df) {
  df %>% select(-circular, -contains("ggraph"))
}


#### ---- File and directory management ---- ####

get_source_data_filename <- function(figure, panel = NULL) {
  final_info <- get_conversion_df() %>%
    filter(make_num == !!figure)

  if (nrow(final_info) != 1) {
    return(NULL)
  } else if (nrow(final_info) > 1) {
    stop(glue("Multiple final numbers found for make figure {figure}."))
  }

  final_num <- final_info$figure_num[[1]]
  is_supp <- final_info$supp[[1]]

  figure_dir <- glue("figure-{str_pad(final_num, 2, pad=\"0\")}")

  if (is_supp) {
    figure_dir <- paste0("supplementary-", figure_dir)
  }

  if (is.null(panel)) {
    filename <- glue("data.tsv")
  } else {
    filename <- glue("panel-{panel}.tsv")
  }

  return(file.path(SOURCE_DATA_BASE_DIR, figure_dir, filename))
}


get_directory <- function(fp) {
  return(str_remove(fp, basename(fp)))
}



#### ---- Data modification ---- ####

round_source_data <- function(x, num_decimals) {
  x %>%
    mutate(across(
      where(is.double),
      scales::label_number(10^-num_decimals)
    ))
}


#### ---- Write data ---- ####

write_source_data <- function(df, filep) {
  dir <- get_directory(filep)
  if (!dir.exists(dir)) {
    dir.create(dir)
  }

  write_tsv(df, filep)
}


save_figure_source_data <- function(x, figure, panel = NULL, num_decimals = 3) {
  if ("ggplot" %in% class(x)) {
    x <- pull_original_plot_data(x)
  } else if (!"tbl_df" %in% class(x)) {
    stop(paste(
      "Cannot save source data for `x`:",
      paste(class(x), collapse = ", ")
    ))
  }

  if (!is.na(num_decimals)) {
    x <- round_source_data(x, num_decimals)
  }

  fn <- get_source_data_filename(figure = figure, panel = panel)
  if (is.null(fn)) {
    warning(glue("Build figure # {figure} is not a final figures (ignoring)"))
    return()
  }
  write_source_data(x, fn)
  invisible(x)
}


#### ---- Archivin ---- ####

archive_source_data <- function() {
  message(glue("{symbol$tick} Archiving Source Data"))
  archive_name <- ""
  archive_path <- here::here("paper", archive_name)
  cmd <- glue("cd paper && tar -czvf source-data.tar.gz source-data && cd ..")
  system(cmd)
}


#### ---- Cleaning ---- ####

clear_source_data <- function() {
  message(glue("{symbol$warning} Clearing cached Source Data directory."))
  if (dir.exists(SOURCE_DATA_BASE_DIR)) {
    unlink(SOURCE_DATA_BASE_DIR, recursive = TRUE)
  }
  dir.create(SOURCE_DATA_BASE_DIR)
}
