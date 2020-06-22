# List of genes with synthetic lethal interactions with KRAS.
# comparisons:
#   1. syn. let. with KRAS (mut vs. WT)
#   2. syn. let. with G12D
#   3. syn. let. with G13D
#   4. syn. let. with G12V


GRAPHS_DIR <- "90_30_synlet-for-shikha"
reset_graph_directory(GRAPHS_DIR)

TABLES_DIR <- GRAPHS_DIR
reset_table_directory(TABLES_DIR)


set.seed(0)


synlet_data <- model1_tib %>%
    filter(cancer == "COAD") %>%
    select(hugo_symbol, data)


gene_is_synlet_possible <- function(df) {
    sum(df$gene_effect < df$cutoff_between_non_essential) >= 2
}

genes_with_synthetic_level_possiblities <- synlet_data %>%
    mutate(syn_let = map_lgl(data, gene_is_synlet_possible)) %>%
    filter(syn_let) %>%
    u_pull(hugo_symbol)

synlet_data %<>%
    filter(hugo_symbol %in% genes_with_synthetic_level_possiblities)

coad_dep_map_ids <- synlet_data %>%
    unnest(data) %>%
    pull(dep_map_id) %>%
    unlist() %>%
    unique()


# Write out table of KRAS mutations in CCLE for Shikha.
ccle_cell_lines %>%
    filter(!is.na(cancer)) %>%
    left_join(ccle_kras_muts, by = "dep_map_id") %>%
    select(dep_map_id, ccle_name, stripped_cell_line_name,
           cancer, lineage, lineage_subtype, lineage_sub_subtype,
           lineage_molecular_subtype, disease, disease_subtype,
           primary_or_metastasis,
           sex, age, culture_medium,
           kras_allele = allele,
           kras_codon = codon,
           kras_cn = copy_number) %>%
    distinct() %>%
    mutate(kras_allele = ifelse(is.na(kras_allele), "WT", kras_allele),
           kras_codon = ifelse(is.na(kras_codon), "WT", kras_codon),
           kras_copynumber = ifelse(is.na(kras_cn), 2, kras_cn)) %>%
    filter(cancer == "COAD") %>%
    write_tsv(table_path(TABLES_DIR, "CCLE-cell-line-information.tsv"))


#### ---- Basic statistical testing ---- ####


# Test one allele vs. all others.
test_allele_vs_other <- function(df, test_allele) {
    df %<>% mutate(kras_allele = allele == !!test_allele)
    t.test(gene_effect ~ kras_allele, data = df)
}


# Get p-value from the results of a t-test.
get_ttest_pval <- function(ttest_res) {
    ttest_res$p.value
}


# Get the mean differences between the groups in the t-test.
get_mean_diff <- function(ttest_res) {
    vals <- unname(ttest_res$estimate)
    vals[2] - vals[1]
}


synlet_results <- synlet_data %>%
    mutate(
        kras_vs_wt = purrr::map(data, test_allele_vs_other,
                                test_allele = "WT"),
        kras_vs_wt_pval = purrr::map_dbl(kras_vs_wt, get_ttest_pval),
        g12d_vs_rest = purrr::map(data, test_allele_vs_other,
                                  test_allele = "G12D"),
        g12d_vs_rest_pval = purrr::map_dbl(g12d_vs_rest, get_ttest_pval),
        g13d_vs_rest = purrr::map(data, test_allele_vs_other,
                                  test_allele = "G13D"),
        g13d_vs_rest_pval = purrr::map_dbl(g13d_vs_rest, get_ttest_pval),
        g12v_vs_rest = purrr::map(data, test_allele_vs_other,
                                  test_allele = "G12V"),
        g12v_vs_rest_pval = purrr::map_dbl(g12v_vs_rest, get_ttest_pval),
    )


# GENES TO USE FOR TESTING COMPLEX MODELS
TEST_GENES <- synlet_results %>%
    filter(g12d_vs_rest_pval < 0.01) %>%
    arrange(g12d_vs_rest_pval) %>%
    slice(1:25) %>%
    select(hugo_symbol, g12d_vs_rest_pval) %T>%
    print() %>%
    pull(hugo_symbol)


#### ---- Model 1. Modeling with comutation interactions ---- ####

# Comutation partners for each allele.
coad_comutations <- genetic_interaction_df %>%
    filter(cancer == "COAD") %>%
    select(allele, hugo_symbol, p_val) %>%
    arrange(allele, p_val, hugo_symbol)
coad_comutations %>%
    count(allele)


# CCLE mutations for cell lines of COAD.
coad_ccle_mutations <- ccle_mutations %>%
    filter(dep_map_id %in% coad_dep_map_ids) %>%
    filter(variant_classification %in% !!coding_mut_var_classes)


# For any numeric columns in a data frame `df`, replace NA's with zero.
replace_numeric_NAs <- function(df) {
    mutate_if(df, is.numeric, replace_na_zero)
}


# Get the comutating genes for a KRAS allele. (memoized)
get_comutating_genes <- function(allele) {
    genes <- coad_comutations %>%
        filter(allele == !!allele) %>%
        pull(hugo_symbol) %>%
        unlist()
}
get_comutating_genes <- memoise::memoise(get_comutating_genes)


# Create a wide tibble of mutations for each cell line in the genes that
# comutate with `allele`. (memoized)
get_comutation_model_data <- function(allele) {
    comut_genes <- get_comutating_genes(allele)

    cna_muts <- coad_ccle_cna

    muts <- coad_ccle_mutations %>%
        filter(hugo_symbol %in% comut_genes) %>%
        add_column(is_mut = 1) %>%
        distinct(dep_map_id, hugo_symbol, is_mut) %>%
        pivot_wider(dep_map_id,
                    names_from = hugo_symbol,
                    values_from = is_mut) %>%
        replace_numeric_NAs()

    allele_comut_plot <- plot_path(
        GRAPHS_DIR,
        as.character(glue("{allele}-ccle-comutation-matrix.svg"))
    )
    if (!file.exists(allele_comut_plot)) {
        plt <- coad_ccle_mutations %>%
            filter(hugo_symbol %in% comut_genes) %>%
            distinct(dep_map_id, hugo_symbol) %>%
            add_column(is_mut = 1) %>%
            complete(dep_map_id, hugo_symbol, fill = list(is_mut = 0)) %>%
            mutate(is_mut = factor(is_mut)) %>%
            ggplot(aes(x = dep_map_id, y = hugo_symbol)) +
            geom_tile(aes(fill = is_mut), color = NA) +
            scale_x_discrete(expand = c(0, 0)) +
            scale_y_discrete(expand = c(0, 0)) +
            scale_fill_manual(values = c("1" = "#6498FF",
                                         "0" = "white")) +
            theme_bw(base_size = 8, base_family = "Arial") +
            theme(
                axis.text.x = element_text(angle = 30, hjust = 1)
            )
        ggsave_wrapper(
            plt,
            allele_comut_plot,
            "large"
        )
    }

    return(muts)
}
get_comutation_model_data <- memoise::memoise(get_comutation_model_data)


# Return the data frame with the covariates except for those of comut. genes.
extract_core_modeling_data <- function(df, allele) {
    df %>%
        select(dep_map_id, gene_effect, rna_scaled, allele, is_altered) %>%
        mutate(allele = as.numeric(allele == !!allele))
}


# Make a model matrix of the core covariates.
core_model_matrix <- function(df) {
    model.matrix(~ 1 + gene_effect + rna_scaled + is_altered, data = df)
}


# Remove the comutation covariates that are the same as `allele`.
remove_covariates_identical_to_allele <- function(mm) {
    allele_vals <- mm[, "allele"]
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
        filter(col_name != "allele") %>%
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
    mm <- model.matrix( ~ -1 + allele * ., data = df) %>%
        remove_covariates_identical_to_allele() %>%
        merge_identical_comutation_covariates()
}


# Tune an fit a `glmnet()` elastic net model using 'caret'. The parameter
# grid is restricted to favor LASSO over Ridge. (memoized)
tune_and_fit_glmnet_elastic <- function(mm) {
    tune_grid <- expand.grid(alpha = seq(0.5, 1, 0.1),
                             lambda = seq(0.0001, 0.5, length = 10))

    message("Tuning elastic net...")
    # Tune the elastic net using 'caret'.
    elastic <- train(
        gene_effect ~ .,
        data = as.data.frame(mm),
        method = "glmnet",
        trControl = trainControl("boot", number = 25),
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


# Model 1.
#   model the gene effect by the RNA expression of the gene, the KRAS allele,
#   mutations to any comutating genes with the allele, and the interactions
#   between the KRAS allele and the comutating gene. The model is fit using
#   using elastic net to find the best alpha and lambda.
model_dependency_with_comutation_1 <- function(df, allele) {
    set.seed(0)

    dat <- extract_core_modeling_data(df, allele) %>%
        left_join(get_comutation_model_data(allele), by = "dep_map_id") %>%
        replace_numeric_NAs() %>%
        select(-dep_map_id)

    # Set up model matrix
    mm1 <- core_model_matrix(dat)
    mm2 <- dat %>%
        select(-c(gene_effect, rna_scaled, is_altered)) %>%
        comutation_model_matrix()
    mm <- cbind(mm1, mm2)

    too_many_zeros <- apply(mm, 2, function(x) { sum(x != 0) < 2 })
    mm <- mm[, !too_many_zeros]

    fit <- tune_and_fit_glmnet_elastic(mm)

    return(list(
        data = dat,
        caret_tune = fit$caret_tune,
        elastic_model = fit$elastic_model,
        alpha = fit$caret_tune$bestTune$alpha,
        lambda = fit$caret_tune$bestTune$lambda,
        allele = allele
    ))

}


# A plot of the coefficient estimates.
coef_plot <- function(mdl) {
    broom::tidy(mdl) %>%
        mutate(term = fct_reorder(term, estimate)) %>%
        ggplot(aes(x = estimate, y = term)) +
        geom_vline(xintercept = 0, lty = 2, size = 0.5, color = "grey25") +
        geom_point() +
        labs(x = "estimate", y = NULL)
}


# Plot the results of model 1.
model_dependency_with_comutation_1_plot <- function(hugo_symbol, fit_results) {

    r3 <- function(x) { round(x, 3) }
    model_caption <- glue(
        "alpha: {r3(fit_results$alpha)}; lambda: {r3(fit_results$lambda)}"
    )

    # coefficient plot
    p_coefs <- coef_plot(fit_results$elastic_model) +
        labs(title = NULL)

    # RNA expression
    rna_pal <- c("grey20", "grey50")
    names(rna_pal) <- c(fit_results$allele, "other")

    p_rna_gene_eff <- fit_results$data %>%
        mutate(
            allele = ifelse(allele == 1, fit_results$allele, "other"),
            is_altered = ifelse(is_altered == 1, "mut", "not. mut")
        ) %>%
        ggplot(aes(x = rna_scaled, y = gene_effect)) +
            geom_point(aes(color = allele, shape = is_altered), size = 2) +
            geom_smooth(method = "lm", formula = "y ~ x") +
            scale_color_manual(values = rna_pal) +
            scale_shape_manual(values = c(17, 16)) +
            labs(x = "RNA expression (scaled)",
                 y = "CERES gene effect",
                 color = "KRAS mut.",
                 shape = glue("{hugo_symbol} mut."))


    # Table of coefficients
    em_coefs_tbl <- broom::tidy(fit_results$elastic_model) %>%
        select(-step, -lambda, -dev.ratio) %>%
        mutate_if(is.numeric, function(x) { round(x, 4) })

    # Box-plots of comutating genes
    comut_terms <- em_coefs_tbl %>%
        filter(!term %in% c("(Intercept)", "allele", "rna_scaled")) %>%
        pull(term)

    comut_terms <- c("allele", comut_terms)

    get_comutation_term_data <- function(df, term) {
        term_split <- unlist(str_split(term, ":|,"))

        if (length(term_split) == 1) {
            mod_df <- df %>%
                select(gene_effect, allele, grp = tidyselect::matches(term))

            if (term == "allele") {
                mod_df$allele <- mod_df$grp
            }
        } else {
            row_sumed <- df %>%
                select(tidyselect::all_of(term_split)) %>%
                apply(1, function(x) { all(x == 1) }) %>%
                as.numeric()
            mod_df <- df %>%
                select(gene_effect, allele) %>%
                mutate(grp = row_sumed)
        }
        return(mod_df)
    }

    facet_nrow <- length(comut_terms) %/% 6 + 1

    comut_terms_plot <- tibble(comutation_term = comut_terms) %>%
        mutate(data = map(comutation_term,
                          get_comutation_term_data,
                          df = fit_results$data)) %>%
        unnest(data) %>%
        mutate(
            allele = ifelse(allele == 1, fit_results$allele, "other"),
            grp = ifelse(grp == 0, "no", "yes")
        ) %>%
        ggplot(aes(x = grp, y = gene_effect)) +
        facet_wrap(~ comutation_term,
                   nrow = facet_nrow,
                   scales = "free") +
        geom_boxplot(fill = NA, color = "grey50", outlier.shape = NA) +
        geom_jitter(aes(color = allele), alpha = 0.7, width = 0.2, size = 0.7) +
        scale_color_manual(values = rna_pal) +
        labs(x = "in comutation group",
             y = "CERES gene effect",
             color = "KRAS mut.")

    patch <- (p_coefs | p_rna_gene_eff ) / comut_terms_plot +
        plot_layout(heights = c(2, facet_nrow)) +
        plot_annotation(title = hugo_symbol, caption = model_caption) &
        theme_bw(base_size = 8, base_family = "Arial")

    ggsave_wrapper(patch,
                   plot_path(GRAPHS_DIR, glue("model1_{hugo_symbol}.svg")),
                   "large")
}



sample_model1 <- synlet_data %>%
    filter(hugo_symbol %in% !!TEST_GENES) %>%
    mutate(
        fit1 = map(data, model_dependency_with_comutation_1, allele = "G12D")
    )

sample_model1 %>%
    mutate(temp_col = map2(
        hugo_symbol, fit1,
        model_dependency_with_comutation_1_plot
    ))



#### ---- Model 2. Including gene effect of comutating genes ---- ####

get_comutation_geneeffect_model_data <- function(allele) {
    comut_genes <- get_comutating_genes(allele)

    muts <- synlet_data %>%
        filter(hugo_symbol %in% comut_genes) %>%
        unnest(data) %>%
        pivot_wider(dep_map_id,
                    names_from = hugo_symbol,
                    values_from = gene_effect) %>%
        replace_numeric_NAs()
    return(muts)
}


model_dependency_with_comutation_2 <- function(df, allele) {
    set.seed(0)

    dat <- extract_core_modeling_data(df, allele) %>%
        left_join(get_comutation_geneeffect_model_data(allele),
                  by = "dep_map_id") %>%
        replace_numeric_NAs() %>%
        select(-dep_map_id)

    # Set up model matrix
    mm <- model.matrix(~ ., data = dat)
    fit <- tune_and_fit_glmnet_elastic(mm)

    return(list(
        data = dat,
        caret_tune = fit$caret_tune,
        elastic_model = fit$elastic_model,
        alpha = fit$caret_tune$bestTune$alpha,
        lambda = fit$caret_tune$bestTune$lambda,
        allele = allele
    ))
}



model_dependency_with_comutation_2_plot <- function(hugo_symbol, fit_results) {
    # coefficient plot
    p_coefs <- coefplot::coefplot(fit_results$elastic_model) +
        theme_bw(base_size = 8, base_family = "Arial") +
        labs(title = hugo_symbol)

    ggsave_wrapper(
        p_coefs,
        plot_path(GRAPHS_DIR, glue("model2_{hugo_symbol}_coefs.svg")),
        size = "medium"
    )

    top_coefs <- broom::tidy(fit_results$elastic_model) %>%
        filter(!str_detect(term, "Intercept|scaled|\\:")) %>%
        top_n(5, wt = abs(estimate)) %>%
        pull(term)

    scatmat_data <- fit_results$data %>%
        select(-rna_scaled) %>%
        mutate(
            allele = ifelse(allele == 1, fit_results$allele, "other"),
            is_altered = ifelse(is_altered == 1, "mut", "not. mut")
        ) %>%
        select(allele, is_altered, gene_effect,
               tidyselect::all_of(top_coefs)) %>%
        dplyr::rename(!!hugo_symbol := gene_effect)

    if (ncol(scatmat_data) > 3) {
        scatmat <- GGally::ggscatmat(scatmat_data, color = "allele") +
            theme_bw(base_size = 9, base_family = "Arial") +
            theme(
                strip.background = element_blank()
            )
        ggsave_wrapper(
            scatmat,
            plot_path(GRAPHS_DIR, glue("model2_{hugo_symbol}_scatmat.svg")),
            "large"
        )
    }
}


sample_model2 <- sample_model1 %>%
    mutate(
        fit2 = map(data, model_dependency_with_comutation_2, allele = "G12D")
    )

sample_model2 %>%
    mutate(temp_col = map2(
        hugo_symbol, fit2,
        model_dependency_with_comutation_2_plot
    ))




#### ---- Box plots ---- ####

# Order the genes by the mean differences of the test.
order_by_mean_diffs <- function(genes, test_res) {
    mean_diffs <- purrr::map_dbl(test_res, get_mean_diff)
    idx <- order(mean_diffs)
    factor(genes, levels = unique(genes[idx]))
}


# Box-plots of `test_allele` vs. the rest.
boxplot_allele_vs_other <- function(df, test_allele) {
    set.seed(0)

    df <- df %>% mutate(
        kras_allele = ifelse(allele == !!test_allele, !!test_allele, "other"
    ))

    pal <- c(short_allele_pal, "other" = "grey75")

    p <- ggplot(df, aes(x = kras_allele, y = gene_effect)) +
        facet_wrap(~ hugo_symbol, scales = "free_y") +
        geom_boxplot(
            aes(color = kras_allele, fill = kras_allele),
            alpha = 0.2,
            outlier.shape = NA,
            width = 0.8
        ) +
        geom_jitter(
            aes(color = allele),
            width = 0.25,
            size = 0.7,
            alpha = 1.0
        ) +
        geom_hline(yintercept = 0, size = 0.5, linetype = 2) +
        scale_color_manual(values = pal) +
        scale_fill_manual(values = pal, guide = FALSE) +
        theme_classic(base_size = 7, base_family = "arial") +
        theme(
            axis.title.x = element_blank(),
            strip.background = element_blank()
        ) +
        labs(y = "depletion effect")
    return(p)
}


boxplot_synlet_results <- function(pval_col, ttest_res_col, test_allele, save_name) {
    pval_col <- rlang::enquo(pval_col)
    ttest_res_col <- rlang::enquo(ttest_res_col)

    synlet_results %>%
        filter(!!pval_col < 0.01) %>%
        mutate(
            hugo_symbol = order_by_mean_diffs(hugo_symbol, !!ttest_res_col)
        ) %>%
        select(hugo_symbol, data) %>%
        unnest(data) %>%
        boxplot_allele_vs_other(test_allele = test_allele) %T>%
        ggsave_wrapper(
            plot_path(GRAPHS_DIR, save_name),
            width = 12, height = 12
        ) %>%
        ggsave_wrapper(
            paste0(file_sans_ext(plot_path(GRAPHS_DIR, save_name)), ".pdf"),
            width = 12, height = 12,
            device = cairo_pdf
        )
}

boxplot_synlet_results(kras_vs_wt_pval, kras_vs_wt,
                       "WT", "mut-vs-other_boxplots.svg")
boxplot_synlet_results(g12d_vs_rest_pval, g12d_vs_rest,
                       "G12D", "G12D-vs-other_boxplots.svg")
boxplot_synlet_results(g13d_vs_rest_pval, g13d_vs_rest,
                       "G13D", "G13D-vs-other_boxplots.svg")
boxplot_synlet_results(g12v_vs_rest_pval, g12v_vs_rest,
                       "G12V", "G12V-vs-other_boxplots.svg")


#### ---- Heatmap ---- ####

make_wide_dataframe <- function(df, scale_rows = FALSE) {
    if (scale_rows) {
        df %<>%
            group_by(hugo_symbol) %>%
            mutate(
                gene_effect = scale(gene_effect)[, 1],
                gene_effect = minmax(gene_effect, -2, 2)
            )
    }
    df %>%
        pivot_wider(id_cols = hugo_symbol,
                    names_from = dep_map_id,
                    values_from = gene_effect) %>%
        as.data.frame() %>%
        column_to_rownames("hugo_symbol")
}

make_row_annotation_dataframe <- function(df, pval_cut) {
    df %>%
        select(hugo_symbol, kras_vs_wt_pval, g12d_vs_rest_pval,
               g13d_vs_rest_pval, g12v_vs_rest_pval) %>%
        unique() %>%
        mutate(
            mut = as.character(kras_vs_wt_pval < pval_cut),
            G12D = as.character(g12d_vs_rest_pval < pval_cut),
            G13D = as.character(g13d_vs_rest_pval < pval_cut),
            G12V = as.character(g12v_vs_rest_pval < pval_cut),
        ) %>%
        select(hugo_symbol, mut, G12D, G12V, G13D) %>%
        as.data.frame() %>%
        column_to_rownames("hugo_symbol")
}


make_col_annotation_dataframe <- function(df) {
    df %>%
        select(dep_map_id, allele) %>%
        unique() %>%
        as.data.frame() %>%
        column_to_rownames("dep_map_id")
}

hm_pal <- colorRampPalette(
    rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu"))
)(11)


clustered_heatmap <- function(df, pval_cut = 0.01, save_name,
                              width, height, fontsize_row = 7,
                              cutree_rows = 0, cutree_cols = 0) {
    wide_df <- make_wide_dataframe(df, scale_rows = TRUE)
    row_anno <- make_row_annotation_dataframe(df, pval_cut)
    col_anno <- make_col_annotation_dataframe(df)

    allele_pal <- c(
        short_allele_pal[names(short_allele_pal) %in% unique(col_anno$allele)],
        "mut" = "black"
    )

    phmat <- pheatmap::pheatmap(
        wide_df,
        annotation_row = row_anno,
        annotation_col = col_anno,
        clustering_distance_rows = "manhattan",
        clustering_distance_cols = "euclidean",
        clustering_method = "ward.D2",
        cutree_rows = cutree_rows,
        cutree_cols = cutree_cols,
        treeheight_row = 60,
        treeheight_col = 20,
        annotation_colors = list(
            mut = c("FALSE" = "white", "TRUE" = "black"),
            G12D = c("FALSE" = "white", "TRUE" = "black"),
            G13D = c("FALSE" = "white", "TRUE" = "black"),
            G12V = c("FALSE" = "white", "TRUE" = "black"),
            allele = allele_pal
        ),
        border_color = NA,
        fontsize = 7,
        fontsize_row = fontsize_row,
        color = hm_pal,
        silent = TRUE
    )

    svg_save_name <- paste0(file_sans_ext(save_name), ".svg")
    pdf_save_name <- paste0(file_sans_ext(save_name), ".pdf")
    save_pheatmap_svg(phmat, plot_path(GRAPHS_DIR, svg_save_name),
                      width = width, height = height)
    save_pheatmap_pdf(phmat, plot_path(GRAPHS_DIR, pdf_save_name),
                      width = width, height = height)
}

{
    # All results on one heatmap.
    synlet_results %>%
        filter(kras_vs_wt_pval < 0.01 |
               g12d_vs_rest_pval < 0.01 |
               g13d_vs_rest_pval < 0.01 |
               g12v_vs_rest_pval < 0.01) %>%
        select(hugo_symbol, data,
               kras_vs_wt_pval, g12d_vs_rest_pval,
               g13d_vs_rest_pval, g12v_vs_rest_pval) %>%
        unnest(data) %>%
        clustered_heatmap(
            save_name = "all-results_heatmap",
            width = 8, height = 15, fontsize_row = 3,
            cutree_rows = 6, cutree_cols = 4
        )

    # KRAS mut vs. WT results.
    synlet_results %>%
        filter(kras_vs_wt_pval < 0.01) %>%
        select(hugo_symbol, data,
               kras_vs_wt_pval, g12d_vs_rest_pval,
               g13d_vs_rest_pval, g12v_vs_rest_pval) %>%
        unnest(data) %>%
        clustered_heatmap(save_name = "mut-results_heatmap",
            width = 6, height = 8, fontsize_row = 7,
            cutree_rows = 2, cutree_cols = 2
        )

    # G12D results.
    synlet_results %>%
        filter(g12d_vs_rest_pval < 0.01) %>%
        select(hugo_symbol, data,
               kras_vs_wt_pval, g12d_vs_rest_pval,
               g13d_vs_rest_pval, g12v_vs_rest_pval) %>%
        unnest(data) %>%
        clustered_heatmap(
            save_name = "g12d-results_heatmap",
            width = 5, height = 10, fontsize_row = 5,
            cutree_rows = 2, cutree_cols = 2
        )

    # G13D results.
    synlet_results %>%
        filter(g13d_vs_rest_pval < 0.01) %>%
        select(hugo_symbol, data,
               kras_vs_wt_pval, g12d_vs_rest_pval,
               g13d_vs_rest_pval, g12v_vs_rest_pval) %>%
        unnest(data) %>%
        clustered_heatmap(
            save_name = "g13d-results_heatmap",
            width = 5, height = 9,
            cutree_rows = 2, cutree_cols = 2
        )

    # G12V results.
    synlet_results %>%
        filter(g12v_vs_rest_pval < 0.01) %>%
        select(hugo_symbol, data,
               kras_vs_wt_pval, g12d_vs_rest_pval,
               g13d_vs_rest_pval, g12v_vs_rest_pval) %>%
        unnest(data) %>%
        clustered_heatmap(
            save_name = "g12v-results_heatmap",
            width = 5, height = 9, fontsize_row = 5,
            cutree_rows = 2, cutree_cols = 2
        )
}


#### ---- Write tables ---- ####

SPREADSHEET <- table_path(TABLES_DIR, "kras-synlet-results.xlsx")

synlet_results_filtered <- synlet_results %>%
    filter(kras_vs_wt_pval < 0.01 |
           g12d_vs_rest_pval < 0.01 |
           g13d_vs_rest_pval < 0.01 |
           g12v_vs_rest_pval < 0.01) %>%
    select(hugo_symbol, data,
           kras_vs_wt_pval, g12d_vs_rest_pval,
           g13d_vs_rest_pval, g12v_vs_rest_pval) %>%
    unnest(data)

synlet_results_filtered %>%
    pivot_wider(
        id_cols = c(hugo_symbol,
                    kras_vs_wt_pval, g12d_vs_rest_pval,
                    g13d_vs_rest_pval, g12v_vs_rest_pval),
        names_from = dep_map_id,
        values_from = gene_effect
    ) %>%
    xlsx::write.xlsx(
        file = SPREADSHEET, sheetName = "data - DepMapID",
        append = TRUE,
    )

synlet_results_filtered %>%
    mutate(allele = paste(allele, dep_map_id, sep = "-")) %>%
    pivot_wider(
        id_cols = c(hugo_symbol,
                    kras_vs_wt_pval, g12d_vs_rest_pval,
                    g13d_vs_rest_pval, g12v_vs_rest_pval),
        names_from = allele,
        values_from = gene_effect
    ) %>%
    xlsx::write.xlsx(
        file = SPREADSHEET, sheetName = "data - allele",
        append = TRUE,
    )

cell_line_names_map <- cell_lines %>%
    select(dep_map_id, stripped_cell_line_name) %>%
    unique() %>%
    dplyr::rename(ccle_name = stripped_cell_line_name) %T>%
    xlsx::write.xlsx(
        file = SPREADSHEET, sheetName = "cell line names",
        append = TRUE,
    )

synlet_results_filtered %>%
    left_join(cell_line_names_map, by = "dep_map_id") %>%
    mutate(cell_line_name = ifelse(
        is.na(ccle_name), dep_map_id, ccle_name
    )) %>%
    pivot_wider(
        id_cols = c(hugo_symbol,
                    kras_vs_wt_pval, g12d_vs_rest_pval,
                    g13d_vs_rest_pval, g12v_vs_rest_pval),
        names_from = cell_line_name,
        values_from = gene_effect
    ) %>%
    xlsx::write.xlsx(
        file = SPREADSHEET, sheetName = "data - cell line names",
        append = TRUE,
    )
