
message("Begninning script.")

.libPaths(c(.libPaths(), "/home/jc604/R-4.0/library"))

message("Loading libraries.")

library(tidygraph)
library(wext)

# Create a bipartite graph from a TSV of samples and their mutations
#     | samples | mutations |


# WRAPPER: make the sample-gene bipartite graph
make_bipartite_graph <- function(data_file, save_file) {
  message("Beginning `make_bipartite_graph()`.")
  message(glue::glue("'data_file': {data_file}"))
  message(glue::glue("'save_file': {save_file}"))

  message("Reading in data.")

  dat <- data.table::fread(data_file)

  message("Making bipartite graph.")

  bgr <- wext::make_sample_gene_bipartite(
    s = unlist(dat[, 1]),
    g = unlist(dat[, 2])
  )

  message("Saving graph.")

  saveRDS(bgr, save_file)
}



# run by Snakemake
make_bipartite_graph(
  data_file = snakemake@input[["in_file"]],
  save_file = snakemake@output[["save_name"]]
)

message("Finished script.")
