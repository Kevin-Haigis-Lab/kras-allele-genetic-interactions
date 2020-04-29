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



#### ---- Statistical testing ---- ####


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



#### ---- Modeling with comutation interactions ---- ####

# Comutation partners for each allele.
coad_comutations <- genetic_interaction_df %>%
    filter(cancer == "COAD" & genetic_interaction == "comutation") %>%
    select(allele, hugo_symbol, p_val) %>%
    arrange(allele, p_val, hugo_symbol)
coad_comutations %>%
    count(allele)

# CCLE mutations for cell lines of COAD.
coad_ccle_mutations <- ccle_mutations %>%
    filter(dep_map_id %in% coad_dep_map_ids) %>%
    filter(variant_classification %in% !!coding_mut_var_classes)


# Create a wide tibble of mutations for each cell line in the genes that
# comutate with `allele`. (memoized)
get_comutation_model_data <- function(allele) {
    comut_genes <- coad_comutations %>%
        filter(allele == !!allele) %>%
        pull(hugo_symbol) %>%
        unlist()

    muts <- coad_ccle_mutations %>%
        filter(hugo_symbol %in% comut_genes) %>%
        add_column(is_mut = 1) %>%
        distinct(dep_map_id, hugo_symbol, is_mut) %>%
        pivot_wider(dep_map_id, names_from = hugo_symbol, values_from = is_mut) %>%
        mutate_if(is.numeric, replace_na_zero)
}
get_comutation_model_data <- memoise::memoise(get_comutation_model_data)


# Model 1.
#   model the gene effect by the RNA expression of the gene, the KRAS allele,
#   mutations to any comutating genes with the allele, and the interactions
#   between the KRAS allele and the comutating gene. The model is fit using
#   using elastic net to find the best alpha and lambda.
model_dependency_with_comutation_1 <- function(df, allele) {
    dat <- df %>%
        select(dep_map_id, gene_effect, rna_scaled, allele, is_altered) %>%
        mutate(allele = as.numeric(allele == !!allele)) %>%
        left_join(get_comutation_model_data(allele), by = "dep_map_id") %>%
        mutate_if(is.numeric, replace_na_zero) %>%
        select(-dep_map_id)

    # Set up model matrix
    mm1 <- model.matrix(~ 1 + gene_effect + rna_scaled + is_altered, data = dat)
    dat_mod <- dat %>% select(-c(gene_effect, rna_scaled, is_altered))
    mm2 <- model.matrix( ~ -1 + allele * ., data = dat_mod)
    mm <- cbind(mm1, mm2)
    too_many_zeros <- apply(mm, 2, function(x) { sum(x != 0) < 2 })
    mm <- mm[, !too_many_zeros]

    message("Tuning elastic net...")
    # Tune the elastic net using 'caret'.
    elastic <- train(
        gene_effect ~ .,
        data = as.data.frame(mm),
        method = "glmnet",
        trControl = trainControl("cv", number = 5),
        tuneLength = 10
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
        data = dat,
        caret_tune = elastic,
        elastic_model = fit,
        alpha = elastic$bestTune$alpha,
        lambda = elastic$bestTune$lambda,
        allele = allele
    ))

}


model_dependency_with_comutation_1_plot <- function(hugo_symbol, fit_results) {
    # coefficient plot
    p_coefs <- coefplot(fit_results$elastic_model) +
        theme_bw(base_size = 8, base_family = "Arial") +
        labs(title = hugo_symbol)

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
            geom_smooth(method = "lm") +
            scale_color_manual(values = rna_pal) +
            scale_shape_manual(values = c(17, 16)) +
            theme_bw(base_size = 8, base_family = "Arial") +
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
        filter(!str_detect(term, ":")) %>%
        pull(term)

    comut_terms <- c("allele", comut_terms)

    get_comutation_term_data <- function(df, term) {
        mod_df <- df %>%
            select(gene_effect, allele, grp = tidyselect::matches(term))
        if (term == "allele") {
            mod_df$allele <- mod_df$grp
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
        theme_bw(base_size = 8, base_family = "Arial") +
            labs(x = "in comutation group",
                 y = "CERES gene effect",
                 color = "KRAS mut.")

    patch <- (p_coefs | p_rna_gene_eff ) / comut_terms_plot +
        plot_layout(heights = c(2, facet_nrow))

    ggsave_wrapper(patch,
                   plot_path(GRAPHS_DIR, glue("{hugo_symbol}_model1.svg")),
                   "large")
}



TEST_GENES <- synlet_results %>%
    arrange(g12d_vs_rest_pval) %>%
    slice(1:10) %>%
    select(hugo_symbol, g12d_vs_rest_pval) %T>%
    print() %>%
    pull(hugo_symbol)

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




sample_fit <- sample_model1$fit1[[2]]$elastic_model
sample_hugo <- sample_model1$hugo_symbol[[2]]
sample_data <- sample_model1$data[[2]]

library(coefplot)
x <- coefplot(sample_fit) +
    theme_bw(base_size = 8,
             base_family = "Arial") +
    labs(title = sample_hugo)
ggsave_wrapper(x, plot_path(GRAPHS_DIR, "test_coefplot.svg"), "medium")


# ggfortify
x <- autoplot(sample_fit, xvar = "lambda")
ggsave_wrapper(x, plot_path(GRAPHS_DIR, "test_autoplot.svg"), "medium")


# broom

# GGally
x <- GGally::ggcoef(sample_fit) +
    theme_bw(base_size = 8, base_family = "Arial") +
    labs(title = sample_hugo)
ggsave_wrapper(x, plot_path(GRAPHS_DIR, "test_ggally_coef.svg"), "medium")


sample_fit
x <- GGally::ggpairs(sample_data)



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
