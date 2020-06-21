# Exploration of the coefficients of all of the models.


GRAPHS_DIR <- "40_65_overview-plot-synlet-comut"
reset_graph_directory(GRAPHS_DIR)

# library(stats)
# library(glue)
# library(conflicted)
# library(ggfortify)
# library(tidygraph)
# library(jhcutils)
# library(magrittr)
# library(ggpubr)
# library(ggraph)
# library(ggtext)
# library(patchwork)
# library(ggplot2)
# library(broom)
# library(tidyverse)
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
# options(dplyr.summarise.inform = FALSE)
#
# load("cache/synlet_comut_model_res.RData")
#
# # synlet_comut_model_res %<>%
# #     filter(
# #         (cancer == "COAD" & allele == "G12D" & hugo_symbol == "STARD9") |
# #             (cancer == "PAAD" & allele == "G12D" & hugo_symbol == "EEF1E1") |
# #             (cancer == "COAD" & allele == "G12D" & hugo_symbol == "SRSF5") |
# #             (cancer == "PAAD" & allele == "G12R" & hugo_symbol == "KIAA1257") |
# #             (cancer == "PAAD" & allele == "G12D" & hugo_symbol == "FKBP1A")
# #     )
#
# load("cache/genetic_interaction_gr.RData")
# load("cache/genetic_interaction_df.RData")
#
#
# saveFigRds <- function(...) { invisible(NULL) }
# reset_graph_directory <- saveFigRds
#
# source("lib/ggplot2_helpers.R")
# source("lib/helpers.R")


# Get the coefficients from a `fit` object and return the results as a tibble.
coef_to_tibble <- function(x) {
    coef(x$elastic_model) %>%
        as.matrix() %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        as_tibble() %>%
        set_names(c("coef", "coef_value")) %>%
        mutate(coef = str_remove_all(coef, "`"))
}

# Coefficients that will often be ignored.
ignore_coefs <- c("(Intercept)", "kras_allele",
                  "rna_expression_std", "is_mutatedTRUE")

# A function to remove these coefficients.
rm_ignore_coefs <- function(df) {
    df %>%
        filter(!coef %in% !!ignore_coefs)
}

# A data frame of the coefficient values.
synlet_comut_coef_data <- synlet_comut_model_res %>%
    mutate(fit_coef = map(fit, coef_to_tibble)) %>%
    select(cancer, allele, hugo_symbol, fit_coef) %>%
    unnest(fit_coef) %>%
    filter(abs(coef_value) > 0.001) %>%
    left_join(
        genetic_interaction_df %>%
            select(hugo_symbol, allele, cancer, genetic_interaction, p_val),
        by = c("coef" = "hugo_symbol", "allele", "cancer")
    )


coefficient_volcano_plot <- synlet_comut_coef_data %>%
    rm_ignore_coefs() %>%
    ggplot(aes(x = coef_value, y = -log10(p_val))) +
    facet_wrap(~ cancer, nrow = 2, scales = "free") +
    geom_point(aes(color = genetic_interaction), size = 2, alpha = 0.5) +
    scale_color_manual(values = comut_mutex_pal) +
    theme_bw(base_size = 7, base_family = "Arial") +
    labs(x = "coefficient value",
         y = "p-value of comutation interaction",
         color = "genetic interaction")
ggsave_wrapper(
    coefficient_volcano_plot,
    plot_path(GRAPHS_DIR, "coefficient_volcano_plot.svg"),
    "small"
)

coefficient_density_plot <- synlet_comut_coef_data %>%
    rm_ignore_coefs() %>%
    ggplot(aes(x = coef_value)) +
    facet_wrap(~ cancer, nrow = 2, scales = "free") +
    geom_density(aes(color = genetic_interaction, fill = genetic_interaction),
                 size = 1, alpha = 0.2) +
    scale_color_manual(values = comut_mutex_pal) +
    scale_fill_manual(values = comut_mutex_pal) +
    theme_bw(base_size = 7, base_family = "Arial") +
    labs(x = "coefficient value",
         y = "density",
         color = "genetic interaction", fill = "genetic interaction")
ggsave_wrapper(
    coefficient_density_plot,
    plot_path(GRAPHS_DIR, "coefficient_density_plot.svg"),
    "small"
)

coefficient_mean_volcano_plot <- synlet_comut_coef_data %>%
    rm_ignore_coefs() %>%
    group_by(cancer, coef, genetic_interaction) %>%
    summarise(mean_coef_value = mean(coef_value),
              mean_p_val = mean(p_val)) %>%
    ungroup() %>%
    ggplot(aes(x = mean_coef_value, y = -log10(mean_p_val))) +
    facet_wrap(~ cancer, nrow = 2, scales = "free") +
    geom_point(aes(color = genetic_interaction), size = 2, alpha = 0.5) +
    scale_color_manual(values = comut_mutex_pal) +
    theme_bw(base_size = 7, base_family = "Arial") +
    labs(x = "coefficient value",
         y = "p-value of comutation interaction",
         color = "genetic interaction")
ggsave_wrapper(
    coefficient_mean_volcano_plot,
    plot_path(GRAPHS_DIR, "coefficient_mean_volcano_plot.svg"),
    "small"
)


coefficient_box_plot <- synlet_comut_coef_data %>%
    rm_ignore_coefs() %>%
    mutate(coef = fct_reorder(coef, coef_value, .fun = mean)) %>%
    ggplot(aes(x = coef, y = coef_value)) +
    facet_wrap(~ cancer, nrow = 2, scales = "free") +
    geom_boxplot(aes(color = genetic_interaction, fill = genetic_interaction),
                 alpha = 0.2) +
    scale_color_manual(values = comut_mutex_pal) +
    scale_fill_manual(values = comut_mutex_pal) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "coefficient",
         y = "coefficient value",
         color = "genetic interaction",
         fill = "genetic interaction")
ggsave_wrapper(
    coefficient_box_plot,
    plot_path(GRAPHS_DIR, "coefficient_box_plot.svg"),
    "wide"
)

# Are the coeffients in the same direction (positive vs. negative)
# as the KRAS allele.
is_same_direction_as_allele <- function(coef, val) {
    if ("kras_allele" %in% coef) {
        kras_allele_val <- val[coef == "kras_allele"]
        allele_direction <- kras_allele_val / abs(kras_allele_val)
        all_direction <- val / abs(val)
        return(all_direction == allele_direction)
    } else {
        return(NA)
    }
}

number_same_direction_allele_bar_plot <- synlet_comut_coef_data %>%
    filter(!(coef %in% ignore_coefs) | coef == "kras_allele") %>%
    group_by(cancer, allele, hugo_symbol) %>%
    mutate(same_dir_as_allele = is_same_direction_as_allele(coef, coef_value)) %>%
    ungroup() %>%
    # group_by(cancer, hugo_symbol, genetic_interaction) %>%
    # summarise(same_dir_as_allele = as.logical(median(same_dir_as_allele))) %>%
    # ungroup() %>%
    count(cancer, genetic_interaction, same_dir_as_allele) %>%
    filter(!is.na(genetic_interaction) & !is.na(same_dir_as_allele)) %>%
    ggplot(aes(x = genetic_interaction, y = n)) +
    facet_wrap(~ cancer, nrow = 1, scales = "free") +
    geom_col(aes(fill = same_dir_as_allele, color = same_dir_as_allele),
             alpha = 0.5, position = "dodge", size = 1) +
    scale_color_manual(values = c("tomato", "dodgerblue")) +
    scale_fill_manual(values = c("tomato", "dodgerblue")) +
    scale_y_continuous(expand = expansion(add = c(0, 2))) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(legend.position = "top") +
    labs(x = "genetic interaction with coefficient",
         y = "count",
         fill = "coefficient is in the same direction as the KRAS allele",
         color = "coefficient is in the same direction as the KRAS allele")
ggsave_wrapper(
    number_same_direction_allele_bar_plot,
    plot_path(GRAPHS_DIR, "number_same_direction_allele_bar_plot.svg"),
    "small"
)
