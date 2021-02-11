
ssrc <- function(n1, n2) {
  all_files <- list.files("src", pattern = "R$", full.names = TRUE)
  all_files <- unlist(all_files)

  n1_pad <- str_pad(as.character(n1), width = 2, side = "left", pad = "0")
  n2_pad <- str_pad(as.character(n2), width = 2, side = "left", pad = "0")
  prefix <- as.character(glue("^{n1_pad}_{n2_pad}"))

  f <- all_files[str_detect(basename(all_files), prefix)]
  if (length(f) == 1) {
    message(glue("Running '{f}'"))
    source(f)
  } else if (length(f) == 0) {
    message(glue("No 'src' file found for '{n1_pad}_{n2_pad}...'"))
  } else {
    message("Multiple files found:")
    message(paste(f, collapse = "\n"))
  }
}
