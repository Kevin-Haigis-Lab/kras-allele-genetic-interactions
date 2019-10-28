# Discriminant analysis to separate the cell lines in the DepMap data.

#### ---- Linear Discriminant Analysis (LDA) ---- ####


caclulate_proportions <- function(lda_model, cancer) {
    (lda_model$svd^2) / sum(lda_model$svd^2)
}

# LASSO-penalized logistic regression on WT vs mutant KRAS
lda_wrapper <- function(df) {
    set.seed(0)

    data <- df %>%
        select(allele, dep_map_id, hugo_symbol, gene_effect) %>%
        unique() %>%
        group_by(hugo_symbol) %>%
        mutate(gene_effect = unlist(scale(gene_effect)[, 1])) %>%
        ungroup() %>%
        pivot_wider(id_cols = c(allele, dep_map_id),
                    names_from = hugo_symbol,
                    values_from = gene_effect) %>%
        select(-dep_map_id)

    data[, -1] <- as.data.frame(data[, -1]) %>%
        DMwR::knnImputation(k = 5, scale = FALSE) %>%
        as_tibble()

    fit_lda <- MASS::lda(allele ~ ., data = data)
    proportions <- caclulate_proportions(fit_lda) * 100 %>% round(1)

    p <- predict(object = fit_lda, newdata = data)$x %>%
        as_tibble() %>%
        mutate(allele = data$allele) %>%
        ggplot(aes(x = LD1, y = LD2)) +
        geom_point(aes(color = allele)) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5)
        ) +
    labs(
        x = glue::glue("LD1 ({proportions[[1]]}%)"),
        y = glue::glue("LD1 ({proportions[[2]]}%)"),
        title = "LDA of DepMap data",
        color = "allele"
    )

    save_path <- plot_path("30_20_discriminant-analysis-syn-let",
                           glue("{cancer}_lda.svg"))
    ggsave_wrapper(p, save_path, medium)

    invisible(fit_lda)
}

# Run a LASSO-penalized logistic regression on WT vs. mutant KRAS
#   for each cancer.
model_data %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(lasso_log_model = map2(data, cancer, lda_wrapper))


# Finished the analysis, but it doesn't show anything interesting.
# Documented on 10/28/2019

# TODO: finish up MVP and comment in notebook
