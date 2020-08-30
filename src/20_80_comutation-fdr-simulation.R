# Simulation for the FDR of the comutation analysis

GRAPHS_DIR <- "20_80_comutation-fdr-simulation"
reset_graph_directory(GRAPHS_DIR)

theme_set(theme_minimal(base_family = "Arial"))


mutation_rate_for_probability <- function(p, omega, tau, alpha, beta) {
    omega * qexp(p, tau) + (1-omega) * invgamma::qinvgamma(p, alpha, beta)
}

crc_mutation_rate_for_probability <- function(p) {
    mutation_rate_for_probability(p, omega = 0.22, tau = 1.52, alpha = 2.23, beta = 0.37)
}


sample_mutation_rate <- function(n, omega, tau, alpha, beta) {
    omega * rexp(n, tau) + (1-omega) * invgamma::rinvgamma(n, alpha, beta)
}


sample_crc_mutation_rate <- function(n) {
    sample_mutation_rate(n, omega = 0.22, tau = 1.52, alpha = 2.23, beta = 0.37)
}


#### ---- Simulation parameters ---- ####

num_simulations <- 1000
tumor_sample_seq <- seq(100, 3000, 100)
# crc_mut_rate_range <- crc_mutation_rate_for_probability(c(0.05, 0.90))
# g1_rate_seq <- seq(crc_mut_rate_range[[1]],
#                    crc_mut_rate_range[[2]],
#                    length.out = 10)
# g2_rate_seq <- seq(0.1, 0.6, 0.1)

broad_mutation_rate_seq <- seq(0.0, 0.8, 0.05)

simulation_grid <- expand.grid(simulation_idx = seq(1, num_simulations),
                               num_tumor_samples = tumor_sample_seq,
                               g1_rate = broad_mutation_rate_seq,
                               g2_rate = broad_mutation_rate_seq) %>%
    as_tibble()
simulation_grid
#> # A tibble: 8,670,000 x 4
#>    simulation_idx num_tumor_samples g1_rate g2_rate
#>             <int>             <dbl>   <dbl>   <dbl>
#>  1              1               100       0       0
#>  2              2               100       0       0
#>  3              3               100       0       0
#>  4              4               100       0       0
#>  5              5               100       0       0
#>  6              6               100       0       0
#>  7              7               100       0       0
#>  8              8               100       0       0
#>  9              9               100       0       0
#> 10             10               100       0       0
#> # … with 8,669,990 more rows


#### ---- Run simulation ---- ####

source(file.path("src", "20_83_comutation-fdr-simulation-shared.R"))

RESET_SIMULATION <- TRUE
if (RESET_SIMULATION) {
    message("Resetting comutation FDR simulation intermediates.")
    unlink(temp_input_dir, recursive = TRUE)
    unlink(temp_output_dir, recursive = TRUE)
}


dir.create(temp_input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(temp_output_dir, showWarnings = FALSE, recursive = TRUE)

submit_script <- file.path("src", "20_81_comutation-fdr-simulation-submit.sh")

write_simulation_job_arrays <- function(job_array_idx, data) {
    input_fname <- simulation_input_file_name(job_array_idx)
    write_tsv(data, input_fname)
}

submit_simulation_job_arrays <- function(job_array_idx, data) {
    if (file.exists(simulation_output_file_name(job_array_idx))) { return() }

    write_simulation_job_arrays(job_array_idx, data)
    cmd <- glue("sbatch {submit_script} {job_array_idx}")
    system(cmd, intern = TRUE)
}

submitted_jobs <- simulation_grid %>%
    mutate(job_array_idx = row_number() %% 1000) %>%
    group_by(job_array_idx) %>%
    nest() %>%
    pmap(submit_simulation_job_arrays) %>%
    unlist()

submitted_jobs[!is.null(submitted_jobs)]


#### ---- Read in simulation results ---- ####

cache("comutation_fdr_simulation_results", {
    simulation_results_files <- list.files(temp_output_dir,
                                           full.names = TRUE,
                                           pattern = "qs$")
    comutation_fdr_simulation_results <- rep(list(NA),
                                             length(simulation_results_files))
    for (i in seq(1, length(simulation_results_files))) {
        f <- simulation_results_files[[i]]
        comutation_fdr_simulation_results[[i]] <- qs::qread(f)
    }
    comutation_fdr_simulation_results <- bind_rows(comutation_fdr_simulation_results)
    return(comutation_fdr_simulation_results)
})

# Check that all simulations were completed and loaded.
stopifnot(nrow(comutation_fdr_simulation_results) == nrow(simulation_grid))



#### ---- Analyze simulation results ---- ####

# Would the simulation pass the project's standard for comutation?
project_significance_check <- function(mut_table, p_value) {

    ## THE FOLLOWING PARAMETERS WERE COPIED FROM "20_35_rc-fisher-comparison.R"
    # p-value thresholds for mutual exclusivity and comutation
    p_val_cut_mutex <- 0.01
    p_val_cut_comut <- 0.01

    # Thresholds for the mutation frequency of other gene for mut. ex. and comut.
    mutfreq_mutex <- 0.02
    mutfreq_comut <- 0.01

    # Threshold for the frequency of comutation between the allele and other gene.
    comutfreq_comut <- 0.10
    ############################################################################

    if (p_value >= p_val_cut_comut) { return(FALSE) }
    if (mut_table[2, 2] == 0) { return(FALSE) }

    num_samples_per_cancer <- sum(mut_table)
    mut_freq_1 <- sum(mut_table[2, ]) / num_samples_per_cancer
    mut_freq_2 <- sum(mut_table[, 2]) / num_samples_per_cancer

    if (mut_freq_1 == 0 | mut_freq_2 == 0) { return(FALSE) }

    comut_freq_1 <- mut_table[2, 2] / sum(mut_table[2, ])
    comut_freq_2 <- mut_table[2, 2] / sum(mut_table[, 2])

    check_1 <- (mut_freq_1 > mutfreq_comut | comut_freq_1 > comutfreq_comut)
    check_2 <- (mut_freq_2 > mutfreq_comut | comut_freq_2 > comutfreq_comut)

    return(check_1 & check_2)

}


simulation_res_stats <- comutation_fdr_simulation_results %>%
    mutate(proj_sig = map2_lgl(mutation_table, p_value,
                               project_significance_check)) %>%
    group_by(g1_rate, g2_rate, num_tumor_samples) %>%
    summarise(freq_of_pval_sig = mean(p_value < 0.01),
              freq_of_proj_sig = mean(proj_sig),
              mean_or = mean(estimate),
              median_of = median(estimate)) %>%
    ungroup()


blue_col <- "#4b68d1"
red_col <- "#d14b4b"
yellow_col <- "#e3df94"

pval_plot <- simulation_res_stats %>%
    ggplot(aes(x = g1_rate, y = g2_rate)) +
    facet_wrap(~ num_tumor_samples) +
    geom_tile(aes(fill = freq_of_pval_sig), color = "white") +
    scale_x_continuous(limits = c(min(broad_mutation_rate_seq),
                                  max(broad_mutation_rate_seq)),
                       expand = c(0, 0)) +
    scale_fill_gradient2(low = blue_col, high = red_col,
                         mid = yellow_col, midpoint = 0.2)
ggsave_wrapper(pval_plot,
               plot_path(GRAPHS_DIR, "pvalue-significance_grid.svg"),
               "large")


proj_cutoffs_plot <- simulation_res_stats %>%
    ggplot(aes(x = g1_rate, y = g2_rate)) +
    facet_wrap(~ num_tumor_samples) +
    geom_tile(aes(fill = freq_of_proj_sig), color = "white") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_distiller(type = "seq", palette = "YlGnBu") +
    theme_minimal(base_size=7, base_family = "Arial") +
    labs(x = "gene 1 mutation rate",
         y = "gene 2 mutation rate",
         fill = "FDR")
ggsave_wrapper(proj_cutoffs_plot,
               plot_path(GRAPHS_DIR, "proj-thresholds_grid.svg"),
               "large")


avg_OR_plot <- simulation_res_stats %>%
    ggplot(aes(x = g1_rate, y = g2_rate)) +
    facet_wrap(~ num_tumor_samples) +
    geom_tile(aes(fill = log(mean_or)), color = "white") +
    scale_x_continuous(limits = c(min(broad_mutation_rate_seq),
                                  max(broad_mutation_rate_seq)),
                       expand = c(0, 0)) +
    scale_fill_gradient2(low = blue_col, high = red_col,
                         mid = yellow_col, midpoint = 0)
ggsave_wrapper(avg_OR_plot,
               plot_path(GRAPHS_DIR, "avg-odds-ratio_grid.svg"),
               "large")
