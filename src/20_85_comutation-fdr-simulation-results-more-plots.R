
library(tidyverse)

libs <- "stats, glue, conflicted, assertr, testthat,
glmnet, parallel, caret, ggfortify, tidygraph, jhcutils,
magrittr, ggpubr, ggraph, ggtext, patchwork, ggplot2, broom, tibble, magrittr,
tidyverse"
libs %>% str_split(",") %>% unlist() %>% str_squish() %>% walk(function(l) {
    library(l, character.only = T)
})

library(mustashe)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("setdiff", "dplyr")
conflict_prefer("intersect", "dplyr")
conflict_prefer("cache", "ProjectTemplate")
conflict_prefer("rename", "dplyr")
conflict_prefer("parLapply", "parallel")
conflict_prefer("which", "Matrix")


for (f in list.files("lib", full.names = TRUE, pattern = "R$")) {
    if (f == "lib/globals.R") { next }
    tryCatch({
        source(f)
    }, error = function(e) {
        message(glue("error in {f}: {e}"))
    })

}


stash("comutation_fdr_simulation_results", {
    load("cache/comutation_fdr_simulation_results.RData")

    comutation_fdr_simulation_results <- comutation_fdr_simulation_results %>%
        filter(num_tumor_samples %in% c(100, 500, 1000, 3000))
})




GRAPHS_DIR <- "20_85_comutation-fdr-simulation-results-more-plots"
reset_graph_directory(GRAPHS_DIR)

source(file.path("src", "20_83_comutation-fdr-simulation-shared.R"))

theme_set(theme_bw(base_size = 7, base_family = "Arial"))

demo_g1_rate <- 0.2
demo_g2_rate <- 0.6
demo_num_tumor_samples <- 100

demo_params_results <- comutation_fdr_simulation_results %>%
    filter(num_tumor_samples == demo_num_tumor_samples) %>%
    filter(near(g1_rate, demo_g1_rate) & near(g2_rate, demo_g2_rate))

num_sims <- max(demo_params_results$simulation_idx)


make_rows_cols_idx <- function(df, arrange_col, num_cols) {
    df %>%
        arrange(-{{ arrange_col }}) %>%
        mutate(row_idx = ceiling(row_number() / num_cols)) %>%
        group_by(row_idx) %>%
        arrange({{ arrange_col }}) %>%
        mutate(col_idx = row_number()) %>%
        ungroup()
}


simulation_grid_plot <- function(df, fill_col) {
    df %>%
        ggplot(aes(col_idx, row_idx)) +
        geom_tile(aes(fill = {{ fill_col }}), color = "grey70") +
        scale_x_continuous(expand = c(0, 0), breaks = c(1, seq(0, 100, 20))) +
        scale_y_continuous(expand = c(0, 0), breaks = seq(1, 10)) +
        coord_fixed(ratio = 3) +
        theme(axis.ticks = element_blank(),
              plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5)) +
        labs(x = NULL,
             y = NULL)
}

demo_colorbar = guide_colorbar(barwidth = unit(3, "mm"),
                               barheight = unit(20, "mm"))

p_title <- glue("Results of {num_sims} simulations with lambdas: {demo_g1_rate} & {demo_g2_rate}")
p_subtitle <- glue("Each cell represents the odds ratio from an individual simulation
                   of {demo_num_tumor_samples} tumor samples.")
demo_or_grid_plot <- demo_params_results %>%
    make_rows_cols_idx(estimate, 100) %>%
    simulation_grid_plot(estimate) +
    scale_fill_gradient(low = "white", high = "red", guide = demo_colorbar) +
    labs(fill = "odds ratio",
         title = p_title,
         subtitle = p_subtitle)
ggsave_wrapper(demo_or_grid_plot,
               plot_path(GRAPHS_DIR, "demo_or-grid-plot.svg"),
               "small")


p_subtitle <- glue("Each cell represents the p-value from an individual simulation
                   of {demo_num_tumor_samples} tumor samples.")
demo_pval_grid_plot <- demo_params_results %>%
    make_rows_cols_idx(-p_value, 100) %>%
    simulation_grid_plot(p_value) +
    scale_fill_gradient(low = "blue", high = "white", guide = demo_colorbar) +
    labs(fill = "p-value",
         title = p_title,
         subtitle = p_subtitle)
ggsave_wrapper(demo_pval_grid_plot,
               plot_path(GRAPHS_DIR, "demo_pval-grid-plot.svg"),
               "small")

p_subtitle <- glue("Each cell represents the if an individual simulation of {demo_num_tumor_samples} tumor samples
                   would result in a false positive")
demo_falsepos_grid_plot <- demo_params_results %>%
    make_rows_cols_idx(-p_value, 100) %>%
    mutate(is_proj_sig = map2_lgl(mutation_table, p_value,
                                  project_significance_check)) %>%
    simulation_grid_plot(is_proj_sig) +
    scale_fill_manual(values = c("white", "black"),
                      labels = c(" true negative", "false positive")) +
    theme(legend.key.size = unit(3, "mm")) +
    labs(x = NULL,
         y = NULL,
         fill = NULL,
         title = p_title,
         subtitle = p_subtitle)
ggsave_wrapper(demo_falsepos_grid_plot,
               plot_path(GRAPHS_DIR, "demo_falsepos-grid-plot.svg"),
               "small")
