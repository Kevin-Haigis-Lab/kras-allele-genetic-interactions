
parse_regions <- function(j) {
  map_dfr(j$regions, function(x) {
    tibble(
      name = x$text,
      start = as.numeric(x$start),
      end = as.numeric(x$end),
      misc = x$metadata$description
    )
  })
}

parse_markups <- function(j) {
  map_dfr(j$markups, function(x) {
    tibble(
      name = x$type,
      start = as.numeric(x$start),
      end = as.numeric(x$start)
    )
  })
}

parse_motifs <- function(j) {
  map_dfr(j$motifs, function(x) {
    tibble(
      name = x$type,
      start = as.numeric(x$start),
      end = as.numeric(x$end)
    )
  })
}

parse_pfam_json <- function(j) {
  bind_rows(
    parse_regions(j),
    parse_markups(j),
    parse_motifs(j)
  )
}

expand_regions <- function(df) {
  df %>%
    mutate(codon = map2(start, end, ~ seq(.x, .y))) %>%
    unnest(codon)
}

extract_gene_name_from_pfam_file <- function(f) {
  a <- str_remove(f, "\\.json$") %>%
    str_split_fixed("_", 2) %>%
    unlist()
  return(a[[2]])
}

pfam_json_to_tibble <- function(f) {
  jsonlite::read_json(f) %>%
    parse_pfam_json() %>%
    expand_regions() %>%
    add_column(hugo_symbol = extract_gene_name_from_pfam_file(f))
}

pfam_files <- list.files(
  file.path("data", "pfam"),
  pattern = "json$",
  full.names = TRUE
)

cache("pfam_data", depends = "pfam_files", {
  pfam_data <- map(pfam_files, pfam_json_to_tibble) %>%
    bind_rows() %>%
    group_by(hugo_symbol) %>%
    nest()
  return(pfam_data)
})
