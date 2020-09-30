
# Count the number of lines of code in a file.
line_count_per_file <- function(f) {
  all_lines <- readLines(f, warn = FALSE) %>%
    unlist() %>%
    str_trim()
  all_lines <- all_lines[str_length(all_lines) > 0]
  comment_line_idx <- str_detect(all_lines, "^#")
  all_lines <- all_lines[!comment_line_idx]
  return(length(all_lines))
}


list_files_in_dir <- function(dir) {
  list.files(dir, full.names = TRUE, pattern = "R$|py$|sh$")
}

# The number of lines of code in this project
comutation_proj_lines_of_code <- function() {

  proj_dirs <- c(
      "src",
      "munge",
      "lib",
      "tests",
      file.path("paper", "figures", "make_figure_scripts")
    )

  d <- tibble(dir = proj_dirs) %>%
    mutate(file = map(dir, list_files_in_dir)) %>%
    unnest(file) %>%
    mutate(num_lines = map_dbl(file, line_count_per_file))

  summ_d <- d %>%
    group_by(dir) %>%
    summarise(
      n_files = n_distinct(file),
      min = min(num_lines),
      q25 = quantile(num_lines, 0.25),
      mean = mean(num_lines),
      median = median(num_lines),
      q75 = quantile(num_lines, 0.75),
      max = max(num_lines),
      total = sum(num_lines)
    ) %>%
    ungroup() %>%
    mutate(dir = basename(dir))

  summ_d %>%
    knitr::kable(digits = 0) %>%
    print()
  cat("\n")
  num_files <- nrow(d)
  num_lines <- sum(d$num_lines)
  message(glue("#> There are {num_lines} lines of code in {num_files} files."))
}

comutation_proj_lines_of_code()
