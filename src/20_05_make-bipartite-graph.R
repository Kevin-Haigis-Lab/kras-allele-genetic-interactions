
.libPaths(c(.libPaths(), "/home/jc604/R-3.5.1/library"))

library(tidygraph)
library(wext)

# Create a bipartite graph from a TSV of samples and their mutations
#     | samples | mutations |


# WRAPPER: make the sample-gene bipartite graph
make_bipartite_graph <- function(data_file, save_file) {
    dat <- data.table::fread(data_file)
    bgr <- wext::make_sample_gene_bipartite(
        s = unlist(dat[, 1]),
        g = unlist(dat[, 2])
    )
    saveRDS(bgr, save_file)
}

# run by Snakemake
make_bipartite_graph(
    data_file = snakemake@input[["in_file"]],
    save_file = snakemake@output[["save_name"]]
)
