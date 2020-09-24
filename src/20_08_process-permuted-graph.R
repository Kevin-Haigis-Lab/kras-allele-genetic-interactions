
.libPaths(c(.libPaths(), "/home/jc604/R-4.0/library"))
library(tictoc)
library(parallel)
library(doParallel)
library(foreach)
library(dplyr)
library(wext)

process_permuted_graph <- function(results_tib_path,
                                   perm_grs_paths,
                                   which_test,
                                   save_name,
                                   n_cores) {

  # check for multiple files
  if (length(perm_grs_paths) < 2) {
    stop("less than 2 files passed through `perm_grs_paths")
  }

  # load the correct test function
  which_test <- stringr::str_to_lower(which_test[[1]])
  if (which_test %in% c("exclusivity", "e")) {
    test_func <- wext::calc_mutex_events
  } else if (which_test %in% c("comutation", "c")) {
    test_func <- wext::calc_comut_events
  } else {
    stop(paste(which_test, "is not an option."))
  }

  # read in results data frame
  results_tib <- readRDS(results_tib_path)

  tic("Process permutations")
  cl <- makeCluster(n_cores, type = "FORK")
  registerDoParallel(cl)
  processed_perms <- foreach(
    p = perm_grs_paths,
    .combine = dplyr::bind_rows
  ) %dopar%
    process_perm(perm_gr_path = p, res_tib = results_tib, f = test_func)
  stopCluster(cl)
  toc()

  tic("Compile results")
  rc_results <- processed_perms %>%
    mutate(gene_sets = purrr::map_chr(gene_sets, unlist_gene_sets)) %>%
    group_by(gene_sets) %>%
    summarise(
      t_BM_ge = sum(t_BM_ge),
      t_AM = t_AM[[1]],
      n_perms = n_distinct(perm_file)
    ) %>%
    ungroup() %>%
    dplyr::mutate(p_val = t_BM_ge / n_perms)
  toc()

  # save the result
  saveRDS(rc_results, save_name)
}


process_perm <- function(perm_gr_path, res_tib, f) {
  perm_gr <- readRDS(perm_gr_path)
  new_res_tib <- purrr::pmap(res_tib, update_results_tib, f = f, bgr = perm_gr) %>%
    bind_rows() %>%
    mutate(perm_file = !!perm_gr_path)
  return(new_res_tib)
}


# update the gene sets in the results tracking tibble using the edge-swapped bgr
update_results_tib <- function(gene_sets, t_BM_ge, t_AM, f, bgr) {
  t_BM <- f(gene_sets, bgr)
  if (t_BM >= t_AM) {
    t_BM_ge <- t_BM_ge + 1
  }
  tibble::tibble(
    gene_sets = list(gene_sets),
    t_BM_ge = t_BM_ge,
    t_AM = t_AM
  )
}


# deterministically turn the list of genes `gs` into a single string
unlist_gene_sets <- function(gs) {
  paste0(sort(unlist(gs)), collapse = " - ")
}


#### ---- INTERACTIONS WITH SNAKEMAKE ---- ####

process_permuted_graph(
  results_tib_path = snakemake@input[["results_df"]],
  perm_grs_paths = snakemake@input[["perm_grs"]],
  which_test = snakemake@wildcards[["which_test"]],
  save_name = snakemake@output[["save_name"]],
  n_cores = snakemake@params[["n_cores"]]
)
