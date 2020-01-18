
set.seed(0)

#### ---- Linear Model 1 ---- ####
# The first model is a simple two-step process. First, the depletion effect is
# regressed on the RNA expression of the target gene. If this fails to explain
# the gene effect (this model is not significant), then an ANOVA and pair-wise
# test is performed over the alleles.


info(logger, "Beginning linear model 1.")


# put the alleles in the correct order of turning into factors
#   if there is WT, put that first
#   all of the other alleles are in alphabetical order
get_allele_factor_levels <- function(als) {
    als <- unique(als)
    if (any(als == "WT")) {
        all_but_wt <- sort(als[als != "WT"])
        return(c("WT", all_but_wt))
    } else {
        return(sort(als))
    }
}


# prepare the data for the first linear model
#   put alleles as factors
#   dummy variable for if the gene is mutated
#   scale the RNA expression to mean of 0 and std. dev. of 1
lm1_prepare_data <- function(tib, hugo_symbol) {
    allele_levels <- get_allele_factor_levels(tib$allele)
    mod_tib <- tib %>%
        ungroup() %>%
        mutate(
            allele = factor(allele, levels = allele_levels),
            is_altered = as.numeric(is_altered),
            rna_scaled = scale(rna_expression)[, 1],
            gene_effect_scaled = scale(gene_effect)[, 1]
        )

    if (any(table(mod_tib$dep_map_id) > 1)) {
        stop(glue("Multiple data points for a cell line in {hugo_symbol}."))
    }

    if (all(is.na(mod_tib$rna_scaled))) { mod_tib$rna_scaled <- 0 }
    return(mod_tib)
}


# linear model on RNA expression
lm_on_rna <- function(data, ...) {
    fit <- lm(gene_effect ~ rna_scaled, data = data)
    return(fit)
}


# does the linear model on RNA have a significant p-value? (returns Boolean)
rna_pvalue_is_significant <- function(pval, cutoff = 0.01) {
    if (is.na(pval)) { return(FALSE) }
    if (pval < cutoff) { return(TRUE) }
    return(FALSE)
}


# ANOVA on the alleles
#   only done if the RNA expression of the gene does not explain its effect
anova_wrapper <- function(data, rna_pvalue, scale_gene_effect = FALSE, ...) {
    if (rna_pvalue_is_significant(rna_pvalue)) { return(NA) }

    if (scale_gene_effect) {
        res <- aov(gene_effect_scaled ~ allele, data = data)
    } else {
        res <- aov(gene_effect ~ allele, data = data)
    }
    return(res)
}


# Kruskal-Wallis rank sum test on the alleles
#   only done if the RNA expression of the gene does not explain its effect
kruskal_wrapper <- function(data, rna_pvalue, ...) {
    if (rna_pvalue_is_significant(rna_pvalue)) { return(NA) }

    res <- kruskal.test(gene_effect ~ allele, data = data)
    return(res)
}


# conduct a pair-wise comparison on each pair of alleles
#   only done if the RNA expression of the gene does not explain its effect
pairwise_wrapper <- function(data, rna_pvalue, scale_gene_effect = FALSE, ...) {
    if (rna_pvalue_is_significant(rna_pvalue)) { return(NA) }

    if (scale_gene_effect) {
        res <- pairwise.t.test(data$gene_effect_scaled,
                               data$allele,
                               p.adjust.method = "BH")
    } else {
        res <- pairwise.t.test(data$gene_effect,
                               data$allele,
                               p.adjust.method = "BH")
    }
    return(res)
}

# conduct the first modeling attempt
model1_tib <- model_data %>%
    filter(!cancer %in% c("MM", "SKCM")) %>%
    filter(!(cancer == "LUAD" & allele == "G13D")) %>%
    group_by(cancer, hugo_symbol) %>%
    filter(sum(gene_effect < cutoff_between_non_essential) >= 2) %>%
    nest() %>%
    mutate(data = purrr::map2(data, hugo_symbol, lm1_prepare_data)) %>%
    ungroup() %>%
    mutate(
        rna_lm_fit = map(data, lm_on_rna),
        rna_pvalue = map_dbl(rna_lm_fit, ~ tidy(.x)$p.value[2]),
        allele_aov = map2(data, rna_pvalue, anova_wrapper),
        allele_pairwise = map2(data, rna_pvalue, pairwise_wrapper)
    )

info(logger, "Caching results of model 1.")
cache("model1_tib")


#### ---- Box-plots ---- ####

# directory for save box-plots
GRAPHS_DIR <- plot_path("10_10_linear-modeling-syn-let_boxplots")
reset_graph_directory(GRAPHS_DIR)

ggproto_save_info <- list(
    COAD = list(fig_num = 4, supp = FALSE)
)

save_proto <- function(gg_obj, save_path, cancer) {
    if (cancer %in% names(ggproto_save_info)) {
        save_info <- ggproto_save_info[[cancer]]
        saveRDS(gg_obj,
                get_fig_proto_path(basename(save_path),
                                   figure_num = save_info$fig_num,
                                   supp = save_info$supp))
    }
}

# plot the results of the first analysis
plot_pairwise_test_results <- function(hugo_symbol, cancer, data,
                                       allele_aov, allele_pairwise,
                                       save_proto = FALSE, ...) {
    if (all(is.na(allele_aov)) | all(is.na(allele_pairwise))) { return() }

    if (tidy(allele_aov)$p.value[[1]] >= 0.01) { return() }

    data <- unique(data)

    stat_tib <- compare_means(
        gene_effect ~ allele, data = data,
        method = "t.test", p.adjust.method = "BH"
    ) %>%
        filter(p.adj < 0.05)

    if (nrow(stat_tib) < 1) { return() }

    stat_bar_height <- 0.08
    stat_bar_y_positions <- c(max(data$gene_effect) + stat_bar_height)
    for (i in seq(1, nrow(stat_tib))) {
        stat_bar_y_positions <- c(
            stat_bar_y_positions,
            stat_bar_y_positions[(i - 1)] + stat_bar_height
        )
    }

    stat_tib$y.position <- stat_bar_y_positions

    p <- ggboxplot(
            data,
            x = "allele",
            y = "gene_effect",
            color = "allele",
            add = "jitter"
        ) +
        stat_pvalue_manual(stat_tib, label = "p.adj", family = "Arial") +
        scale_color_manual(values = short_allele_pal) +
        theme(
            text = element_text(family = "arial"),
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_blank(),
            legend.position = "none"
        ) +
        labs(
            title = glue("{hugo_symbol} in {cancer}"),
            caption = "*bars indicate FDR-adjusted p-values < 0.05",
            y = "depletion effect"
        )

    plot_fname <- file.path(GRAPHS_DIR, glue("{cancer}-{hugo_symbol}.svg"))
    ggsave_wrapper(p, plot_fname, size = "small")

    if (save_proto) {
        save_proto(p, plot_fname, cancer)
    }
}


select_gene_boxplots <- tibble::tribble(
    ~ cancer, ~ hugo_symbol,
      "COAD",        "IDH1",
      "COAD",       "KNTC1",
      "COAD",     "PIP5K1A",
      "COAD",       "WDR26",
)

model1_tib %>%
    pwalk(plot_pairwise_test_results) %>%
    right_join(select_gene_boxplots, by = c("cancer", "hugo_symbol")) %>%
    pwalk(plot_pairwise_test_results, save_proto = TRUE)



#### ---- Heatmaps ---- ####

cancer_pheatmap_manager <- list(
    COAD = list(
        col_cuts = 4,
        row_cuts = 5
    ),
    LUAD = list(
        col_cuts = 2,
        row_cuts = 4
    ),
    PAAD = list(
        col_cuts = 3,
        row_cuts = 6
    )
)


prep_pheatmap_df <- function(data, method = c("scale", "normalize")) {
    f <- function(x) { x }
    if (method == "scale") {
        f <- function(x) { scale(x)[, 1] }
    } else if (method == "normalize") {
        f <- function(x) { scales::rescale(x, to = c(-2, 2)) }
    } else {
        stop(glue("method not a possible option: {method}"))
    }

    df <- data %>%
        group_by(hugo_symbol) %>%
        mutate(gene_effect_scaled = f(gene_effect)) %>%
        ungroup() %>%
        select(hugo_symbol, dep_map_id, gene_effect_scaled) %>%
        unique() %>%
        group_by(hugo_symbol, dep_map_id) %>%
        filter(n() == 1) %>%
        ungroup() %>%
        pivot_wider(
            names_from = dep_map_id,
            values_from = gene_effect_scaled) %>%
        column_to_rownames("hugo_symbol")
    invisible(df)
}


merge_wt <- function(df, data) {
    wt_samples <- data %>%
        filter(allele == "WT") %>%
        pull(dep_map_id) %>%
        unique()
    wt_df <- apply(df[, wt_samples], 1, mean, na.rm = TRUE)
    new_df <- df[, !colnames(df) %in% wt_samples]
    new_df$WT <- wt_df
    return(new_df)
}

fig_save_info <- list(
    COAD = list(fig_num = 4, supp = FALSE)
)

save_proto <- function(cancer, ph, save_path) {
    if (cancer %in% names(fig_save_info)) {
        save_info <- fig_save_info[[cancer]]
    } else {
        return(NULL)
    }

    saveRDS(ph, get_fig_proto_path(basename(save_path),
                                   save_info$fig_num,
                                   supp = save_info$supp)
    )
}


plot_cancer_heatmaps <- function(cancer, data, screen,
                                 merge_luad = TRUE,
                                 row_dist_method = "euclidean",
                                 col_dist_method = "euclidean",
                                 row_hclust_method = "complete",
                                 col_hclust_method = "complete") {

    mod_data <- prep_pheatmap_df(data, "normalize")

    if (cancer == "LUAD" & merge_luad) {
        mod_data <- merge_wt(df = mod_data, data = data)
    }

    row_hclust <- hclust(dist(mod_data, method = row_dist_method),
                         method = row_hclust_method)
    col_hclust <- hclust(dist(t(mod_data), method = col_dist_method),
                         method = col_hclust_method)

    col_anno <- data %>%
        select(dep_map_id, allele) %>%
        unique() %>%
        column_to_rownames("dep_map_id")

    row_anno <- cutree(
            row_hclust,
            k = cancer_pheatmap_manager[[cancer]]$row_cuts
        ) %>%
        enframe(value = "cluster") %>%
        mutate(cluster = factor(cluster)) %>%
        column_to_rownames("name")

    if (cancer == "LUAD") {
        wt_anno_df <- data.frame(allele = "WT")
        rownames(wt_anno_df) <- "WT"
        col_anno <- rbind(col_anno, wt_anno_df)
    }

    anno_pal <- list(allele = short_allele_pal[as.character(unique(col_anno$allele))])

    pal <- c(synthetic_lethal_pal["down"], "grey95", synthetic_lethal_pal["up"])
    pal <- colorRampPalette(pal)(7)

    ph <- pheatmap::pheatmap(
        mod_data,
        color = pal,
        cluster_rows = row_hclust,
        cluster_cols = col_hclust,
        annotation_col = col_anno,
        annotation_row = row_anno,
        annotation_colors = anno_pal,
        cutree_rows = cancer_pheatmap_manager[[cancer]]$row_cuts,
        cutree_cols = cancer_pheatmap_manager[[cancer]]$col_cuts,
        treeheight_row = 8,
        treeheight_col = 8,
        fontsize = 5,
        show_rownames = (cancer != "LUAD"),
        silent = TRUE,
        border_color = NA,
        fontfamily = "Arial"
    )

    save_path <- plot_path(
        "10_10_linear-modeling-syn-let_pheatmaps",
        glue("{cancer}_{screen}_{row_dist_method}_{row_hclust_method}_pheatmap.svg")
    )
    save_pheatmap_svg(ph, save_path, width = 7, height = 9)
    save_proto(cancer, ph, save_path)
}


cluster_genes <- function(cancer, data,
                          row_dist_method = "euclidean",
                          row_hclust_method = "complete") {
    mod_data <- prep_pheatmap_df(data, method = "normalize")

    if (cancer == "LUAD") {
        mod_data <- merge_wt(df = mod_data, data = data)
    }

    gene_hclust <- hclust(dist(mod_data, method = row_dist_method),
                          method = row_hclust_method)
    gene_cls <- cutree(
            gene_hclust,
            k = cancer_pheatmap_manager[[cancer]]$row_cuts
        ) %>%
        enframe(name = "hugo_symbol", value = "gene_cls")

    return(gene_cls)
}


# Run `plot_cancer_heatmaps()` with a bunch of distance and clustering methods.
plot_cancer_heatmaps_multiple_methods <- function(cancer, data, screen,
                                                  merge_luad = TRUE,
                                                  methods_df) {
    pwalk(methods_df, plot_cancer_heatmaps,
          cancer = cancer, data = data,
          screen = screen, merge_luad = merge_luad)
}


# Combinations of `dist()` and `hclust()` methods.
dist_methods <- c("euclidean", "manhattan")
hclust_methods <- c("ward.D", "ward.D2", "single", "complete",
                    "average", "mcquitty", "median", "centroid")
methods_tib <- expand.grid(dist_methods, hclust_methods,
                           stringsAsFactors = FALSE) %>%
    as_tibble()
colnames(methods_tib) <- c("row_dist_method", "row_hclust_method")


# make a heatmap for the genes in each cancer
depmap_gene_clusters <- model1_tib %>%
    filter(rna_pvalue > 0.01) %>%
    mutate(
        aov_p_val = purrr::map_dbl(allele_aov, ~tidy(.x)$p.value[[1]])
    ) %>%
    filter(aov_p_val < 0.01) %>%
    select(hugo_symbol, cancer, data) %>%
    unnest(data) %>%
    group_by(cancer) %>%
    nest() %T>%
    purrr::pwalk(plot_cancer_heatmaps_multiple_methods,
             screen = "CRISPR", methods_df = methods_tib) %>%
    mutate(cluster_tib = purrr::map2(cancer, data, cluster_genes,
                                     row_dist_method = "manhattan",
                                     row_hclust_method = "ward.D2")) %>%
    select(-data) %>%
    unnest(cluster_tib) %>%
    ungroup()

cache("depmap_gene_clusters", depends = "model1_tib")





#### ---- RNAi Screen modeling ---- ####


# conduct the first modeling attempt
rnai_model1_tib <- rnai_model_data %>%
    group_by(cancer, allele, hugo_symbol) %>%
    filter(n_distinct(dep_map_id) >= 3) %>%
    group_by(cancer, hugo_symbol) %>%
    filter(n_distinct(allele) > 1) %>%
    nest() %>%
    mutate(data = purrr::map2(data, hugo_symbol, lm1_prepare_data)) %>%
    ungroup() %>%
    mutate(
        rna_lm_fit = map(data, lm_on_rna),
        rna_pvalue = map_dbl(rna_lm_fit, ~ tidy(.x)$p.value[2]),
        allele_aov = map2(data, rna_pvalue, anova_wrapper),
        allele_pairwise = map2(data, rna_pvalue, pairwise_wrapper)
    )

rnai_model1_tib %>%
    filter(rna_pvalue > 0.01) %>%
    mutate(aov_p_val = purrr::map_dbl(allele_aov, ~ tidy(.x)$p.value[[1]])) %>%
    filter(aov_p_val < 0.01) %>%
    select(hugo_symbol, cancer, data) %>%
    unnest(data) %>%
    group_by(cancer) %>%
    nest() %>%
    ungroup() %T>%
    purrr::pwalk(plot_cancer_heatmaps, screen = "RNAi")


info(logger, "Caching results of RNAi model 1.")
cache("rnai_model1_tib")
