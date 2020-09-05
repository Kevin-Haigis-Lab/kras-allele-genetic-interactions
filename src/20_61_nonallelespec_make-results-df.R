
.libPaths(c(.libPaths(), "/home/jc604/R-3.5.1/library"))

library(tictoc)
library(tidygraph)
library(dplyr)
library(wext)

# A modified RC-test from the 'wext' package optimized for this analysis.
mod_rc_test <- function(bipartite_gr,
                        which_test = c("exclusivity", "comutation"),
                        k = 2,
                        seed_genes = c(),
                        min_mut_events = 3,
                        min_times_mut = 2) {
    # choose only one test
    which_test <- stringr::str_to_lower(which_test[[1]])
    if (which_test %in% c("exclusivity", "e")) {
        test_func <- wext::calc_mutex_events
    } else if (which_test %in% c("comutation", "c")) {
        test_func <- wext::calc_comut_events
    } else {
        stop(paste(which_test, "is not an option."))
    }

    # load bipartite graph (real)
    bipartite_gr <- readRDS(bipartite_gr)

    # create the empty results tibble
    tic("created results tibble")
    results_tib <- make_empty_results_tracker(bipartite_gr, k,
                                              which_test = which_test,
                                              seed_genes = seed_genes,
                                              min_times_mut = min_times_mut) %>%
        dplyr::mutate(
            t_AM = purrr::map_dbl(gene_sets, test_func, bgr = bipartite_gr)
        ) %>%
        dplyr::filter(t_AM >= !!min_mut_events)
    toc()

    return(results_tib)
}



#### ---- INTERACTIONS WITH SNAKEMAKE ---- ####

# WRAPPER: run the RC test with the real and permuted bipartite graphs
run_rc_test <- function(real_gr,
                        which_test,
                        min_times_mut,
                        output_name) {

    if (seed_genes == "NULL") seed_genes <- c()
    stopifnot(length(seed_genes) > 0)

    results <- mod_rc_test(bipartite_gr = real_gr,
                           which_test = which_test,
                           k = 2,
                           seed_genes = c(seed_genes),
                           min_mut_events = 2,
                           min_times_mut = min_times_mut)

    saveRDS(results, output_name)
}

# adds an underscore between the RAS and allele
convert_rasallele <- function(rasallele) {
    ras <- stringr::str_extract(rasallele, "KRAS|NRAS")
    allele <- stringr::str_extract(rasallele, "(?<=RAS)[:alnum:]+")
    return(paste0(ras, "_", allele))
}


# run by Snakemake
run_rc_test(
    real_gr = snakemake@input[["real_gr"]],
    which_test = snakemake@wildcards[["which_test"]],
    min_times_mut = as.numeric(snakemake@params[["min_times_mut"]]),
    output_name = snakemake@output[["output_name"]]
)
