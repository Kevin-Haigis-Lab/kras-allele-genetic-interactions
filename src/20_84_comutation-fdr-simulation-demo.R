
# library(tidyverse)
#
# libs <- "stats, glue, conflicted, assertr, testthat,
# glmnet, parallel, caret, ggfortify, tidygraph, jhcutils,
# magrittr, ggpubr, ggraph, ggtext, patchwork, ggplot2, broom, tibble, magrittr,
# tidyverse"
# libs %>% str_split(",") %>% unlist() %>% str_squish() %>% walk(function(l) {
#     library(l, character.only = T)
# })
#
# conflict_prefer("select", "dplyr")
# conflict_prefer("filter", "dplyr")
# conflict_prefer("slice", "dplyr")
# conflict_prefer("setdiff", "dplyr")
# conflict_prefer("intersect", "dplyr")
# conflict_prefer("cache", "ProjectTemplate")
# conflict_prefer("rename", "dplyr")
# conflict_prefer("parLapply", "parallel")
# conflict_prefer("which", "Matrix")
#
#
# for (f in list.files("lib", full.names = TRUE, pattern = "R$")) {
#     tryCatch({
#         source(f)
#     }, error = function(e) {
#         message(glue("error in {f}: {e}"))
#     })
#
# }


GRAPHS_DIR <- "20_84_comutation-fdr-simulation-demo"
reset_graph_directory(GRAPHS_DIR)

source(file.path("src", "20_83_comutation-fdr-simulation-shared.R"))

set.seed(0)

theme_set(
    theme_bw(base_size = 7, base_family = "Arial") %+replace%
        theme(plot.tag = element_text(face = "bold"))
)

lambda_1 <- 0.2
lambda_2 <- 0.6


pois_dist_data <- tibble(lambda = c(lambda_1, lambda_2)) %>%
    mutate(num_events = map(lambda_1, ~ seq(0, 5)),
           prob_density = map2(num_events, lambda, ~ dpois(x = .x, lambda = .y))) %>%
    unnest(c(num_events, prob_density)) %>%
    mutate(lambda_lbl = glue("lambda: {lambda}"))

pois_dist_plots <- pois_dist_data %>%
    ggplot(aes(x = num_events, y = prob_density)) +
    facet_wrap(~ lambda_lbl, nrow = 1) +
    geom_col(fill = "black") +
    scale_x_continuous(breaks = seq(0, 10)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 9, face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    labs(x = "number of mutations",
         y = "probability mass",
         tag = "a")


mut_dist_plots <- pois_dist_data %>%
    mutate(any_events = num_events > 0) %>%
    group_by(lambda, lambda_lbl, any_events) %>%
    summarise(prob_density = sum(prob_density)) %>%
    ungroup() %>%
    ggplot(aes(x = any_events, y = prob_density)) +
    facet_wrap(~ lambda_lbl, nrow = 1) +
    geom_col(aes(fill = any_events), color = "black") +
    scale_x_discrete(labels = c("WT", "mutant")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = c("white", "black"), guide = FALSE) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 9, face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    labs(x = "number of mutations",
         y = "probability mass",
         tag = "b")


N_TS <- 100
simulation_data <- tibble(tumor_sample = seq(1, N_TS)) %>%
    mutate(g1_mut = rpois(N_TS, lambda = lambda_1) > 0,
           g2_mut = rpois(N_TS, lambda = lambda_1) > 0)

comutation_grid_plot <- simulation_data %>%
    mutate(comut = g1_mut & g2_mut) %>%
    pivot_longer(-c(tumor_sample, comut),
                 names_to = "gene",
                 values_to = "is_mut") %>%
    mutate(fill_lbl = case_when(
        comut ~ "comutation",
        is_mut ~ "mutated",
        TRUE ~ "not mutated"
    )) %>%
    mutate(gene_lbl = ifelse(gene == "g1_mut", "gene 1", "gene 2")) %>%
    ggplot(aes(tumor_sample, gene_lbl)) +
    geom_tile(aes(fill = fill_lbl), color = "grey50") +
    scale_x_continuous(expand = c(0, 0), breaks = c(1, seq(0, N_TS, 20))) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = c("blue", "black", "white")) +
    theme(axis.ticks = element_blank(),
          axis.text.y = element_text(face = "bold"),
          legend.position = "bottom",
          legend.key.size = unit(3, "mm"),
          plot.title = element_text(hjust = 0.5)) +
    labs(x = "simulated tumor sample",
         y = NULL,
         title = glue("Mutation events in {N_TS} simulated tumors"),
         fill = NULL,
         tag = "c")


simulated_ct <- my_table(simulation_data$g1_mut, simulation_data$g2_mut)
simulated_fisher <- fisher.test(simulated_ct, alternative = "g")
is_false_positive <- project_significance_check(
    mut_table = simulated_ct,
    p_value = simulated_fisher$p.value
)

simulation_res_tbl <- broom::tidy(simulated_fisher) %>%
    transmute(estimate = round(estimate, 3),
              p.value = round(p.value, 3)) %>%
    rename(`odds ratio` = estimate,
           `p-value` = p.value) %>%
    add_column(`is false positive?` = ifelse(is_false_positive, "yes", "no"))

simulation_res_tbl <- gridExtra::tableGrob(
    simulation_res_tbl,
    rows = NULL,
    theme = gridExtra::ttheme_default(base_size = 8)
)

my_arrow <- arrow(angle = 25, length = unit(2, "mm"), type = "closed")
plt_arrow <- tibble(x = seq(1, 2), y = 1) %>%
    ggplot(aes(x, y)) +
    geom_line(arrow = my_arrow) +
    theme_void()


plt_bent_arrow1 <- tibble(x = c(1, 2, 2), y = c(2, 2, 1)) %>%
    ggplot(aes(x, y)) +
    geom_line(arrow = my_arrow) +
    scale_x_continuous(limits = c(1, 3), expand = c(0, 0)) +
    scale_y_continuous(limits = c(1, 3), expand = c(0, 0)) +
    theme_void()


plt_bent_arrow2 <- tibble(x = c(1, 1, 2), y = c(2, 1, 1)) %>%
    ggplot(aes(x, y)) +
    geom_line(arrow = my_arrow) +
    scale_x_continuous(limits = c(0, 2), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 2), expand = c(0, 0)) +
    theme_void()




simulation_demo_patch <- (
    (
        (plot_spacer() | pois_dist_plots | plt_arrow | mut_dist_plots | plt_bent_arrow1) +
            plot_layout(widths = c(1, 2, 1, 2, 1))
    ) /
    comutation_grid_plot /
    (plt_bent_arrow2 | simulation_res_tbl | plot_spacer())
) +
    plot_layout(heights = c(2, 1, 1))

ggsave_wrapper(simulation_demo_patch,
               plot_path(GRAPHS_DIR, "simulation-demo-patch.svg"),
               "wide")
