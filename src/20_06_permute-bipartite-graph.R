
.libPaths(c(.libPaths(), "/home/jc604/R-4.0/library"))

library(tidygraph)
library(dplyr)
library(wext)

# WRAPPER: run the bipartite edge swapping algorithm
run_bipartite_edge_swap <- function(bgr_file, Q, save_file) {
    bgr <- readRDS(bgr_file)
    perm_gr <- wext::bipartite_edge_swap3(bgr, Q)
    saveRDS(perm_gr, save_file)
}

# run by Snakemake
run_bipartite_edge_swap(
    snakemake@input[["bipartite_gr"]],
    snakemake@params[["Q"]],
    snakemake@output[["save_name"]]
)
