# Explaining the KRAS allele-specific synthetic lethal interactions by
# the comutating genes.

GRAPHS_DIR <- "40_60_synlet-explained-by-comuts"
reset_graph_directory(GRAPHS_DIR)



#### ---- Data preparation subroutines ---- ####

# Get the genes that comutate with an allele in a cancer.
get_comutating_genes <- function(cancer, allele) {
    genetic_interaction_df %>%
        filter(cancer == !!cancer & allele == !!allele) %>%
        pull(hugo_symbol) %>%
        unlist() %>%
        unique()
}


# Get a wide data framge of comutations in the cell lines. (memoised)
get_comutation_dataframe <- function(cancer, allele, cell_lines, min_muts = 0) {
    comut_genes <- get_comutating_genes(cancer, allele)
    ccle_mutations_dmg %>%
        filter(dep_map_id %in% !!cell_lines) %>%
        filter(hugo_symbol %in% !!comut_genes) %>%
        distinct(dep_map_id, hugo_symbol) %>%
        add_column(is_mut = 1) %>%
        add_count(hugo_symbol) %>%
        filter(n >= min_muts) %>%
        select(-n) %>%
        pivot_wider(dep_map_id,
                    names_from = hugo_symbol,
                    values_from = is_mut) %>%
        replace_numeric_NAs()
}
get_comutation_dataframe <- memoise::memoise(get_comutation_dataframe)



#### ---- Model Matrix subroutines ---- ####

# Remove the comutation covariates that are the same as `allele`.
remove_covariates_identical_to_allele <- function(mm) {
    allele_vals <- mm[, "kras_allele"]
    same_as_allele <- apply(mm, 2, function(x) all(x == allele_vals))
    same_as_allele[1] <- FALSE
    return(mm[, !same_as_allele])
}


# Merge all identical comutation covariates into a single covariate.
merge_identical_comutation_covariates <- function(mm) {

    merged_columns <- tibble(
        col_name = colnames(mm),
        values = apply(mm, 2, function(x) {paste0(x, collapse = ",")})
    ) %>%
        filter(col_name != "kras_allele") %>%
        group_by(values) %>%
        summarise(comb_cols = list(col_name),
                  new_col_name = paste0(col_name, collapse = ",")) %>%
        select(new_col_name, comb_cols)


    for (i in seq(1, nrow(merged_columns))) {
        mod_mm <- mm

        rm_cols <- unlist(merged_columns$comb_cols[i])
        rm_cols_idx <- Matrix::which(colnames(mod_mm) %in% rm_cols)
        mod_mm <- mod_mm[, -c(rm_cols_idx)]

        new_col_name <- merged_columns$new_col_name[[i]]
        new_mat <- as.matrix(mm[, rm_cols[[1]]])
        colnames(new_mat) <- new_col_name
        mod_mm <- cbind(mod_mm, new_mat)

        mm <- mod_mm
    }

    return(mm)
}


# Make the model matrix for the comutation covariates.
comutation_model_matrix <- function(df) {
    mm <- model.matrix( ~ -1 + kras_allele * ., data = df) %>%
        remove_covariates_identical_to_allele() %>%
        merge_identical_comutation_covariates()
}


# Make a model matrix of the core covariates.
core_model_matrix <- function(df) {
    model.matrix(~ 1 + gene_effect + rna_expression_std + is_mutated, data = df)
}


# Return the data frame with the covariates except for those of comut. genes.
extract_core_modeling_data <- function(df, allele) {
    df %>%
        select(dep_map_id, gene_effect, rna_expression_std,
               kras_allele, is_mutated) %>%
        mutate(kras_allele = as.numeric(kras_allele == !!allele))
}


#### ---- Fitting Elastic Net ---- ####

# Return a tibble of the results of the best elastic net model.
get_best_glmnet_result <- function(caret_fit) {
    best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
    best_result = caret_fit$results[best, ]
    rownames(best_result) = NULL
    as_tibble(best_result)
}


# Tune an fit a `glmnet()` elastic net model using 'caret'. The parameter
# grid is restricted to favor LASSO over Ridge. (memoized)
tune_and_fit_glmnet_elastic <- function(mm) {
    tune_grid <- expand.grid(alpha = seq(0.75, 1, 0.05),
                             lambda = seq(0.0001, 0.5, length = 20))

    message("Tuning elastic net...")
    # Tune the elastic net using 'caret'.
    elastic <- train(
        gene_effect ~ .,
        data = as.data.frame(mm),
        method = "glmnet",
        trControl = trainControl("boot", number = 30),
        tuneGrid = tune_grid
    )

    # Modify model matrix for `glmnet()`.
    y <- mm[, "gene_effect"]
    mm <- mm[, !(colnames(mm) %in% c("(Intercept)", "gene_effect"))]

    message("Fitting final elastic net...")
    fit <- glmnet(x = mm,
                  y = y,
                  alpha = elastic$bestTune$alpha,
                  lambda = elastic$bestTune$lambda)

    message("Done!")
    return(list(
        caret_tune = elastic,
        elastic_model = fit
    ))
}
tune_and_fit_glmnet_elastic <- memoise::memoise(tune_and_fit_glmnet_elastic)


# Model the dependency score using the standard covariates and vairbales for
# whether a comutating gene is mutated or not.
synlet_with_comutations <- function(cancer, allele, hugo_symbol, data,
                                    min_muts = 3) {
    set.seed(0)

    comutation_df <- get_comutation_dataframe(cancer, allele,
                                              unique(data$dep_map_id),
                                              min_muts = min_muts)

    model_data <- extract_core_modeling_data(data, allele) %>%
        left_join(comutation_df, by = "dep_map_id") %>%
        replace_numeric_NAs() %>%
        select(-dep_map_id)

    mm1 <- core_model_matrix(model_data)
    mm2 <- model_data %>%
        select(-c(gene_effect, rna_expression_std, is_mutated)) %>%
        comutation_model_matrix()
    mm <- cbind(mm1, mm2)

    too_many_zeros <- apply(mm, 2, function(x) { sum(x != 0) < 3 })
    mm <- mm[, !too_many_zeros]

    fit <- tune_and_fit_glmnet_elastic(mm)

    return(list(
        data = model_data,
        caret_tune = fit$caret_tune,
        elastic_model = fit$elastic_model,
        alpha = fit$caret_tune$bestTune$alpha,
        lambda = fit$caret_tune$bestTune$lambda,
        result_diagnostics = get_best_glmnet_result(fit$caret_tune),
        allele = allele
    ))
}


#### ---- Diagnositc and model plotting ---- ####

# A plot of the coefficient estimates.
coef_plot <- function(mdl) {
    broom::tidy(mdl) %>%
        mutate(term = str_replace_all(term, ",", "\n"),
               term = fct_reorder(term, estimate)) %>%
        ggplot(aes(x = estimate, y = term)) +
        geom_vline(xintercept = 0, lty = 2, size = 0.5, color = "grey25") +
        geom_point() +
        labs(x = "estimate", y = NULL)
}


# Plot RNA expression against gene effect.
rna_plot <- function(hugo_symbol, fit_results) {
    p <- fit_results$data %>%
        mutate(
            allele = ifelse(kras_allele == 1, fit_results$allele, "other"),
            allele = factor(allele, levels = c(fit_results$allele, "other")),
            is_mutated = ifelse(is_mutated == 1, "mut", "not. mut"),
            is_mutated = factor(is_mutated, levels = c("mut", "not. mut"))
        ) %>%
        ggplot(aes(x = rna_expression_std, y = gene_effect)) +
        geom_point(aes(color = allele, shape = is_mutated), size = 2) +
        geom_smooth(method = "lm", formula = "y ~ x") +
        scale_color_manual(values = c("grey10", "grey60"),
                           drop = FALSE,
                           guide = guide_legend(order = 5)) +
        scale_shape_manual(values = c(17, 16),
                           drop = FALSE,
                           guide = guide_legend(order = 10)) +
        labs(x = "RNA expression (scaled)",
             y = "CERES gene effect",
             color = "KRAS mut.",
             shape = glue("{hugo_symbol} mut."))
    return(p)
}


# Get the information for whether a cell line was parts of a term or not.
get_comutation_term_data <- function(df, term) {
    term_split <- unlist(str_split(term, ":|,"))

    if (length(term_split) == 1) {
        mod_df <- df %>%
            select(gene_effect, kras_allele, grp = tidyselect::matches(term))

        if (term == "allele") {
            mod_df$kras_allele <- mod_df$grp
        }
    } else {
        row_sumed <- df %>%
            select(tidyselect::all_of(term_split)) %>%
            apply(1, function(x) { all(x == 1) }) %>%
            as.numeric()
        mod_df <- df %>%
            select(gene_effect, kras_allele) %>%
            mutate(grp = row_sumed)
    }
    return(mod_df)
}


# Box-plots for all of the terms in the model.
comut_term_boxplots <- function(fit_results) {
    # Table of coefficients
    em_coefs_tbl <- broom::tidy(fit_results$elastic_model) %>%
        select(-step, -lambda, -dev.ratio) %>%
        mutate_if(is.numeric, function(x) { round(x, 4) })

    comut_terms <- em_coefs_tbl %>%
        filter(
            !term %in% c("(Intercept)", "kras_allele", "rna_expression_std")
        ) %>%
        pull(term)

    comut_terms <- c("allele", comut_terms)

    facet_nrow <- length(comut_terms) %/% 6 + 1

    p <- tibble(comutation_term = comut_terms) %>%
        mutate(data = map(comutation_term,
                          get_comutation_term_data,
                          df = fit_results$data)) %>%
        unnest(data) %>%
        mutate(
            kras_allele = ifelse(kras_allele == 1, fit_results$allele, "other"),
            kras_allele = factor(kras_allele,
                                 levels = c(fit_results$allele, "other")),
            grp = ifelse(grp == 0, "no", "yes")
        ) %>%
        ggplot(aes(x = grp, y = gene_effect)) +
        facet_wrap(~ comutation_term,
                   nrow = facet_nrow,
                   scales = "free") +
        geom_boxplot(fill = NA, color = "grey50", outlier.shape = NA) +
        geom_jitter(aes(color = kras_allele),
                    alpha = 0.7, width = 0.2, size = 0.7) +
        scale_color_manual(values =  c("grey10", "grey60"), drop = FALSE) +
        labs(x = "in comutation group",
             y = "CERES gene effect",
             color = "KRAS mut.")

    return(list(plot = p, num_rows = facet_nrow))
}


# Plot the results of model 1.
synlet_with_comutations_plots <- function(hugo_symbol, fit_results) {

    r3 <- function(x) { round(x, 3) }
    fit_diagnotics <- fit_results$result_diagnostics
    model_caption <- paste(
        glue("alpha: {r3(fit_results$alpha)}"),
        glue("lambda: {r3(fit_results$lambda)}"),
        glue("Rsquared: {r3(fit_diagnotics$Rsquared[[1]])}"),
        glue("Rsquared std. dev.: {r3(fit_diagnotics$RsquaredSD[[1]])}"),
        glue("RMSE: {r3(fit_diagnotics$RMSE[[1]])}"),
        glue("RMSE std. dev.: {r3(fit_diagnotics$RMSESD[[1]])}"),
        glue("MAE: {r3(fit_diagnotics$MAE[[1]])}"),
        glue("MAE std. dev.: {r3(fit_diagnotics$MAESD[[1]])}"),
        sep = "; "
    )

    # coefficient plot
    p_coefs <- coef_plot(fit_results$elastic_model) +
        labs(title = NULL)

    p_rna_gene_eff <- rna_plot(hugo_symbol, fit_results)

    # Box-plots of comutating genes
    comut_terms_plots <- comut_term_boxplots(fit_results)
    facet_nrow <- comut_terms_plots$num_rows
    comut_terms_plots <- comut_terms_plots$plot

    patch <- (p_coefs | p_rna_gene_eff ) / comut_terms_plots +
        plot_layout(heights = c(2, facet_nrow)) +
        plot_annotation(title = hugo_symbol, caption = model_caption) &
        theme_bw(base_size = 8, base_family = "Arial")

    ggsave_wrapper(patch,
                   plot_path(GRAPHS_DIR, glue("synlet-comut_{hugo_symbol}.svg")),
                   "large")

    return(patch)
}



#### ---- Run model ---- ####


cache("synlet_comut_model_res",
      depends = c("depmap_model_workflow_res",
                  "ccle_mutations_dmg",
                  "genetic_interaction_df"),
{
    synlet_comut_model_res <- depmap_model_workflow_res %>%
        filter_depmap_model_workflow_res() %>%
        select(cancer, hugo_symbol, data, ova_pairs) %>%
        unnest(ova_pairs) %>%
        filter(adj_p_value < 0.05 & allele != "WT") %>%
        select(cancer, allele, hugo_symbol, data) %>%
        mutate(fit = pmap(., synlet_with_comutations),
               plt = map2(hugo_symbol, fit, synlet_with_comutations_plots))
    return(synlet_comut_model_res)
})
