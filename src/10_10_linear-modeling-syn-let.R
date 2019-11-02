
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
lm1_prepare_data <- function(tib) {
    allele_levels <- get_allele_factor_levels(tib$allele)
    mod_tib <- tib %>%
        group_by(dep_map_id) %>%
        slice(1) %>%
        ungroup() %>%
        mutate(
            allele = factor(allele, levels = allele_levels),
            is_altered = as.numeric(is_altered),
            rna_scaled = scale(rna_expression)[, 1]
        )
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
anova_wrapper <- function(data, rna_pvalue, ...) {
    if (rna_pvalue_is_significant(rna_pvalue)) { return(NA) }

    res <- aov(gene_effect ~ allele, data = data)
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
pairwise_wrapper <- function(data, rna_pvalue, ...) {
    if (rna_pvalue_is_significant(rna_pvalue)) { return(NA) }

    res <- pairwise.t.test(data$gene_effect,
                           data$allele,
                           p.adjust.method = "BH")
    return(res)
}

# conduct the first modeling attempt
model1_tib <- model_data %>%
    filter(!cancer %in% c("MM", "SKCM")) %>%
    filter(!(cancer == "LUAD" & allele == "G13D")) %>%
    group_by(cancer, hugo_symbol) %>%
    filter(sum(gene_effect < cutoff_between_non_essential) >= 2) %>%
    nest() %>%
    mutate(data = purrr::map(data, lm1_prepare_data)) %>%
    ungroup() %>%
    mutate(
        rna_lm_fit = map(data, lm_on_rna),
        rna_pvalue = map_dbl(rna_lm_fit, ~ tidy(.x)$p.value[2]),
        allele_aov = map2(data, rna_pvalue, anova_wrapper),
        allele_pairwise = map2(data, rna_pvalue, pairwise_wrapper)
    )

info(logger, "Caching results of model 1.")
cache("model1_tib")


# plot the results of the first analysis
plot_pairwise_test_results <- function(hugo_symbol, cancer, data,
                                       allele_aov, allele_pairwise, ...) {
    if (all(is.na(allele_aov)) | all(is.na(allele_pairwise))) { return() }

    if (tidy(allele_aov)$p.value[[1]] >= 0.01) { return() }

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
        stat_pvalue_manual(stat_tib, label = "p.adj") +
        scale_color_manual(values = short_allele_pal) +
        theme(
            text = element_text(family = "arial"),
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_blank(),
            legend.position = "none"
        ) +
        labs(
            title = glue("Cancer: {cancer}, gene: {hugo_symbol}"),
            caption = "*bars indicate FDR-adjusted p-values < 0.05",
            y = "depletion effect"
        )

    plot_fname <- file.path(plot_save_dir, glue("{cancer}-{hugo_symbol}.svg"))
    ggsave_wrapper(p, plot_fname, size = "small")
}


# directory for save images
plot_save_dir <- plot_path("10_10_linear-modeling-syn-let_boxplots")
if (!dir.exists(plot_save_dir)) {
    info(logger, glue("Making directory for saving boxplots: {plot_save_dir}"))
    dir.create(plot_save_dir)
} else {
    all_files <- list.files(plot_save_dir, pattern = "svg", full.names = TRUE)
    file.remove(all_files)
}

model1_tib %>%
    pwalk(plot_pairwise_test_results)



#### ---- Heatmaps ---- ####

cancer_pheatmap_manager <- list(
    COAD = list(
        col_cuts = 4,
        row_cuts = 4
    ),
    LUAD = list(
        col_cuts = 2,
        row_cuts = 4
    ),
    PAAD = list(
        col_cuts = 3,
        row_cuts = 4
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

plot_cancer_heatmaps <- function(cancer, data) {

    df <- prep_pheatmap_df(data, "normalize")

    if (cancer == "LUAD") {
        df <- merge_wt(df = df, data = data)
    }

    row_hclust <- hclust(dist(df, method = "manhattan"))
    col_hclust <- hclust(dist(t(df), method = "manhattan"))

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

    anno_pal <- list(allele = short_allele_pal[unique(col_anno$allele)])

    ph <- pheatmap::pheatmap(
        df,
        cluster_rows = row_hclust,
        cluster_cols = col_hclust,
        annotation_col = col_anno,
        annotation_row = row_anno,
        annotation_colors = anno_pal,
        cutree_rows = cancer_pheatmap_manager[[cancer]]$row_cuts,
        cutree_cols = cancer_pheatmap_manager[[cancer]]$col_cuts,
        treeheight_col = 20,
        show_rownames = (cancer != "LUAD"),
        silent = TRUE
    )

    save_path <- plot_path("10_10_linear-modeling-syn-let_pheatmaps",
                           glue("{cancer}_pheatmap.svg"))
    save_pheatmap_svg(ph, save_path, width = 7, height = 9)
}


cluster_genes <- function(cancer, data) {
    df <- prep_pheatmap_df(data, method = "normalize")

    if (cancer == "LUAD") {
        df <- merge_wt(df = df, data = data)
    }

    gene_hclust <- hclust(dist(df, method = "manhattan"))
    gene_cls <- cutree(
            gene_hclust,
            k = cancer_pheatmap_manager[[cancer]]$row_cuts
        ) %>%
        enframe(name = "hugo_symbol", value = "gene_cls")

    return(gene_cls)
}


# make a heatmap for the genes in each cancer
depmap_gene_clusters <- model1_tib %>%
    filter(rna_pvalue > 0.01) %>%
    mutate(aov_p_val = purrr::map_dbl(allele_aov, ~ tidy(.x)$p.value[[1]])) %>%
    filter(aov_p_val < 0.01) %>%
    select(hugo_symbol, cancer, data) %>%
    unnest(data) %>%
    group_by(cancer) %>%
    nest() %T>%
    purrr::pwalk(plot_cancer_heatmaps) %>%
    mutate(cluster_tib = purrr::map2(cancer, data, cluster_genes)) %>%
    select(-data) %>%
    unnest(cluster_tib) %>%
    ungroup()

cache("depmap_gene_clusters")
