# Discriminant analysis to separate the cell lines in the DepMap data.

#### ---- Linear Discriminant Analysis (LDA) ---- ####


caclulate_proportions <- function(lda_model, cancer) {
  (lda_model$svd^2) / sum(lda_model$svd^2)
}

# LASSO-penalized logistic regression on WT vs mutant KRAS
lda_wrapper <- function(cancer, data, ...) {
  set.seed(0)

  mod_data <- data %>%
    select(allele, dep_map_id, hugo_symbol, gene_effect) %>%
    unique() %>%
    group_by(hugo_symbol) %>%
    mutate(gene_effect = unlist(scale(gene_effect)[, 1])) %>%
    ungroup() %>%
    pivot_wider(
      id_cols = c(allele, dep_map_id),
      names_from = hugo_symbol,
      values_from = gene_effect
    ) %>%
    select(-dep_map_id)

  mod_data[, -1] <- as.data.frame(mod_data[, -1]) %>%
    DMwR::knnImputation(k = 5, scale = FALSE) %>%
    as_tibble()

  fit_lda <- MASS::lda(allele ~ ., data = mod_data)
  proportions <- round(caclulate_proportions(fit_lda) * 100, 1)

  p <- predict(object = fit_lda, newdata = mod_data)$x %>%
    as_tibble() %>%
    mutate(allele = mod_data$allele) %>%
    ggplot(aes(x = LD1, y = LD2)) +
    geom_point(aes(color = allele)) +
    scale_color_manual(values = short_allele_pal) +
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

  save_path <- plot_path(
    "30_20_discriminant-analysis-syn-let",
    glue("{cancer}_lda.svg")
  )
  ggsave_wrapper(p, save_path, "medium")

  invisible(fit_lda)
}

# Run a LASSO-penalized logistic regression on WT vs. mutant KRAS
#   for each cancer.
model_data %>%
  filter(!cancer %in% c("SKCM", "MM")) %>%
  group_by(cancer) %>%
  nest() %>%
  mutate(fit_lda = purrr::map2(cancer, data, lda_wrapper))


# Finished the analysis, but it doesn't show anything interesting.
# Documented on 10/28/2019
