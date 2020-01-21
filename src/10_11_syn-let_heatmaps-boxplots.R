# Box-plots and heatmaps of results from linear modeling of allele-specific
# synthetic lethality.


#### ---- Box-plots ---- ####

# directory for save box-plots
GRAPHS_DIR_BOXES <- "10_11_linear-modeling-syn-let_boxplots"
reset_graph_directory(GRAPHS_DIR_BOXES)

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
            add = "jitter",
            size = ifelse(save_proto, 0.01, 1)
        ) +
        stat_pvalue_manual(stat_tib, label = "p.adj", family = "Arial") +
        scale_color_manual(values = short_allele_pal) +
        theme(
            text = element_text(family = "arial"),
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_blank(),
            legend.position = "none"
        ) +
        labs(y = "depletion effect")

    if (!save_proto) {
        p <- p +
            labs(
                title = glue("{hugo_symbol} in {cancer}"),
                caption = "*bars indicate FDR-adjusted p-values < 0.05"
            )
    }

    plot_fname <- plot_path(GRAPHS_DIR_BOXES, glue("{cancer}-{hugo_symbol}.svg"))
    ggsave_wrapper(p, plot_fname, size = "small")

    if (save_proto) {
        save_proto(p, plot_fname, cancer)
    }
}


# Plot the results of the first analysis.
# This version of the function uses my own box-plot creation function
#   located in "lib/stats-boxplot.R".
plot_pairwise_test_results2 <- function(hugo_symbol, cancer, data,
                                        allele_aov, allele_pairwise,
                                        save_proto = FALSE, replace_svg = FALSE,
                                        ...) {
    if (all(is.na(allele_aov)) | all(is.na(allele_pairwise))) { return() }

    if (tidy(allele_aov)$p.value[[1]] >= 0.01) { return() }

    data <- unique(data) %>%
        mutate(x = fct_drop(factor_alleles(allele)),
               y = gene_effect)


    stat_tib <- make_stats_dataframe(data, auto_filter = TRUE,
                                     method = "t.test", p.adjust.method = "BH")

    if (nrow(stat_tib) < 1) { return() }
    p <- stats_boxplot(data, stat_tib, box_color = allele,
                       up_spacing = 0.06,
                       point_size = 0.1, label_size = 3, bar_size = 0.4) +
        scale_color_manual(values = short_allele_pal) +
        theme_bw(base_size = 7, base_family = "Arial") +
        theme(
            text = element_text(family = "arial"),
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_blank(),
            legend.position = "none"
        ) +
        labs(y = "depletion effect")

    if (replace_svg) {
        plot_fname <- file.path(GRAPHS_DIR_BOXES, glue("{cancer}-{hugo_symbol}.svg"))
        ggsave_wrapper(p, plot_fname, size = "small")
    }

    if (save_proto) {
        save_proto(p, plot_fname, cancer)
    }
}

# Select genes for figures.
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
    pwalk(plot_pairwise_test_results2, save_proto = TRUE)



#### ---- Heatmaps ---- ####

# directory for save box-plots
GRAPHS_DIR_HEAT <- plot_path("10_11_linear-modeling-syn-let_pheatmaps")
reset_graph_directory(GRAPHS_DIR_HEAT)

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
        GRAPHS_DIR_HEAT,
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


# RNAi heatmaps
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
