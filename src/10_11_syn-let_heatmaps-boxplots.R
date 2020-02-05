# Box-plots and heatmaps of results from linear modeling of allele-specific
# synthetic lethality.


#### ---- Box-plots ---- ####

# directory for save box-plots
GRAPHS_DIR_BOXES <- "10_11_linear-modeling-syn-let_boxplots"
reset_graph_directory(GRAPHS_DIR_BOXES)

ggproto_save_info <- list(
    COAD = list(fig_num = 4, supp = FALSE),
    PAAD = list(fig_num = 13, supp = TRUE)
)

save_boxplot_proto <- function(gg_obj, save_path, cancer) {
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
                                       ...) {
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
            size = 1
        ) +
        stat_pvalue_manual(stat_tib, label = "p.adj", family = "Arial") +
        scale_color_manual(values = short_allele_pal) +
        theme(
            text = element_text(family = "Arial"),
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_blank(),
            legend.position = "none"
        ) +
        labs(y = "depletion effect")

    plot_fname <- plot_path(GRAPHS_DIR_BOXES,
                            glue("{cancer}-{hugo_symbol}.svg"))
    ggsave_wrapper(p, plot_fname, size = "small")
}


# Plot the results of the first analysis.
# This version of the function uses my own box-plot creation function
#   located in "lib/stats-boxplot.R".
plot_pairwise_test_results2 <- function(hugo_symbol, cancer, data,
                                        filter_aov = TRUE, filter_stats = TRUE,
                                        allele_aov, allele_pairwise,
                                        save_proto = FALSE,
                                        replace_svg = FALSE,
                                        ...) {
    if (all(is.na(allele_aov)) | all(is.na(allele_pairwise))) return()

    if (filter_aov && tidy(allele_aov)$p.value[[1]] >= 0.01) return()

    data <- unique(data) %>%
        mutate(x = fct_drop(factor_alleles(allele)),
               y = gene_effect)


    stat_tib <- make_stats_dataframe(data, auto_filter = TRUE,
                                     method = "t.test", p.adjust.method = "BH")

    if (filter_stats && nrow(stat_tib) < 1) return()

    p <- stats_boxplot(data, stat_tib, box_color = allele,
                       up_spacing = 0.06,
                       point_size = 0.1, label_size = 3, bar_size = 0.3) +
        scale_color_manual(values = short_allele_pal) +
        theme_bw(base_size = 7, base_family = "Arial") +
        theme(
            text = element_text(family = "Arial"),
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_blank(),
            legend.position = "none"
        ) +
        labs(y = "depletion effect")

    plot_fname <- plot_path(GRAPHS_DIR_BOXES,
                            glue("{cancer}-{hugo_symbol}.svg"))
    if (replace_svg) {
        ggsave_wrapper(p, plot_fname, size = "small")
    }

    if (save_proto) {
        save_boxplot_proto(p, plot_fname, cancer)
    }
}

# Select genes for figures.
select_gene_boxplots <- tibble::tribble(
    ~ cancer, ~ hugo_symbol,
      "COAD",        "IDH1",
      "COAD",       "KNTC1",
      "COAD",     "PIP5K1A",
      "COAD",       "WDR26",
      "PAAD",         "JUN",
      "PAAD",      "CEP350",
      "PAAD",       "EGLN2",
      "PAAD",        "NUMB",
      "PAAD",        "TMED2"
)

model1_tib %>%
    pwalk(plot_pairwise_test_results) %>%
    right_join(select_gene_boxplots, by = c("cancer", "hugo_symbol")) %>%
    pwalk(plot_pairwise_test_results2, save_proto = TRUE)


## Specifically make plot for MAPK8 in PAAD.
model_data %>%
    filter(cancer == "PAAD" & hugo_symbol == "MAPK8") %>%
    group_by(hugo_symbol, cancer) %>%
    nest() %>%
    mutate(allele_aov = TRUE, allele_pairwise = TRUE) %>%
    pmap(plot_pairwise_test_results2,
         filter_aov = FALSE, filter_stats = FALSE,
         save_proto = TRUE, replace_svg = TRUE,)
#####


#### ---- Heatmaps ---- ####

# directory for save box-plots
GRAPHS_DIR_HEAT <- "10_11_linear-modeling-syn-let_pheatmaps"
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
    new_df["WT (avg.)"] <- wt_df
    return(new_df)
}

fig_save_info <- list(
    COAD = list(fig_num = 4, supp = FALSE),
    LUAD = list(fig_num = 4, supp = FALSE),
    PAAD = list(fig_num = 13, supp = TRUE)
)


replace_WT_in_luad <- function(df) {
    df %>% mutate(name = ifelse(name == "WT (avg.)", "WT", name))
}

# Custom annotation legends using `ggplot2::geom_tile()`.
make_annotation_tile <- function(v, grp) {
    df <- tibble::enframe(v) %>% filter(!is.na(name))

    if (grp == "heat") {
        df %<>% mutate(name = as.numeric(name))
    } else {
        df %<>%
            replace_WT_in_luad() %>%
            mutate(name = factor(name, levels = rev(sort(unique(name)))))
    }

    p <- df %>%
        ggplot(aes(x = "", y = name)) +
        geom_tile(aes(fill = value), color = NA) +
        scale_fill_identity() +
        labs(y = grp)
    return(p)
}


# Save a proto for a figure.
save_pheatmap_proto <- function(cancer, ph, save_path,
                                anno_pal = NULL, heat_pal = NULL) {
    if (cancer %in% names(fig_save_info)) {
        save_info <- fig_save_info[[cancer]]
    } else {
        return(NULL)
    }

    saveRDS(ph, get_fig_proto_path(basename(save_path),
                                   save_info$fig_num,
                                   supp = save_info$supp))

    anno_save_path <- function(pal_name) {
        paste0(file_sans_ext(basename(save_path)), "_", pal_name, ".svg")
    }

    if (!is.null(anno_pal)) {
        g_allele <- make_annotation_tile(anno_pal$allele, "allele")
        p <- anno_save_path("allelepal")
        saveRDS(g_allele,
                get_fig_proto_path(p, save_info$fig_num, supp = save_info$supp))
        g_cls <- make_annotation_tile(anno_pal$cluster, "cluster")
        p <- anno_save_path("clusterpal")
        saveRDS(g_cls,
                get_fig_proto_path(p, save_info$fig_num, supp = save_info$supp))
    }

    if (!is.null(heat_pal)) {
        g_heat <- make_annotation_tile(heat_pal, "heat")
        p <- anno_save_path("heatpal")
        saveRDS(g_heat,
                get_fig_proto_path(p, save_info$fig_num, supp = save_info$supp))
    }
}

# Mappings from default cluster assignments to order shown in pheatmap.
cluster_number_map <- list(
    COAD = tibble::tribble(
        ~default_cluster, ~cluster,
        1, 2,
        2, 5,
        3, 1,
        4, 4,
        5, 3,
    ),
    LUAD = tibble::tribble(
        ~default_cluster, ~cluster,
        1, 1,
        2, 4,
        3, 3,
        4, 2,
    ),
    PAAD = tibble::tribble(
        ~default_cluster, ~cluster,
        1, 6,
        2, 2,
        3, 4,
        4, 1,
        5, 5,
        6, 3,
    )
)


# Update the clusters of the final data frame to those manually assigned in
# `cluster_number_map`.
update_clusters <- function(tib) {
    mapping_df <- tibble::enframe(cluster_number_map, name = "cancer") %>%
        unnest(value)

    tib %>%
        dplyr::rename(default_cluster = gene_cls) %>%
        left_join(mapping_df, by = c("cancer", "default_cluster")) %>%
        select(-default_cluster) %>%
        dplyr::rename(gene_cls = cluster)
}


# Get color palette for clusters. It is derived from viridis.
cluster_color_pal <- function(n_vals) {
    cols <- viridis::viridis_pal()(n_vals)
    names(cols) <- seq(1, n_vals)
    return(cols)
}


apriori_genes <- unique(unlist(c(
    kegg_geneset_df$hugo_symbol,
    cosmic_cgc_df$hugo_symbol,
    kras_interactors_bioid_df$hugo_symbol
)))



# Custom row names of only a priori genesets.
labels_for_apriori_genes <- function(df) {
    ifelse(rownames(df) %in% apriori_genes, rownames(df), "")
}


# Plot pretty heatmaps for a cancer.
plot_cancer_heatmaps <- function(cancer, data, screen,
                                 merge_luad = TRUE,
                                 row_dist_method = "euclidean",
                                 col_dist_method = "euclidean",
                                 row_hclust_method = "complete",
                                 col_hclust_method = "complete",
                                 save_proto = TRUE) {
    set.seed(0)

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
        enframe(value = "default_cluster") %>%
        left_join(cluster_number_map[[cancer]], by = "default_cluster") %>%
        mutate(cluster = factor(cluster, levels = sort(unique(cluster)))) %>%
        select(-default_cluster) %>%
        column_to_rownames("name")

    if (cancer == "LUAD") {
        wt_anno_df <- data.frame(allele = "WT (avg.)")
        rownames(wt_anno_df) <- "WT (avg.)"
        col_anno <- rbind(col_anno, wt_anno_df)
    }

    anno_pal <- list(
        allele = short_allele_pal[as.character(sort(unique(col_anno$allele)))],
        cluster = cluster_color_pal(max(as.numeric(row_anno$cluster)))
    )

    if (cancer == "LUAD") {
        idx <- names(anno_pal$allele) == "WT"
        names(anno_pal$allele)[idx] <- "WT (avg.)"
    }

    pal <- c(synthetic_lethal_pal["down"], "grey95", synthetic_lethal_pal["up"])
    pal <- colorRampPalette(pal)(7)

    if (cancer == "LUAD") {
        labels_row <- labels_for_apriori_genes(mod_data)
    } else {
        labels_row <- NULL
    }

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
        labels_row = labels_row,
        silent = TRUE,
        border_color = NA,
        fontfamily = "Arial"
    )

    save_path <- plot_path(
        GRAPHS_DIR_HEAT,
        glue("{cancer}_{screen}_{row_dist_method}_{row_hclust_method}_pheatmap.svg")
    )
    save_pheatmap_svg(ph, save_path, width = 4, height = 5)

    names(pal) <- seq(min(mod_data, na.rm = TRUE),
                      max(mod_data, na.rm = TRUE),
                      length.out = length(pal))
    if (save_proto) {
        save_pheatmap_proto(cancer, ph, save_path,
                            anno_pal = anno_pal, heat_pal = pal)
    }
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
    ungroup() %>%
    update_clusters()

cache("depmap_gene_clusters", depends = "model1_tib")


# Print out the number of genes per cancer.
depmap_gene_clusters %>%
    group_by(cancer) %>%
    summarise(num_genes = n_distinct(hugo_symbol)) %>%
    ungroup()


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
    purrr::pwalk(plot_cancer_heatmaps, screen = "RNAi", save_proto = FALSE)
