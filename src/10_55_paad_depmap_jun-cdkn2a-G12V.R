# Investigating the relationship between JUN and KRAS G12V in PAAD.

GRAPHS_DIR <- "90_25_paad_depmap_jun-cdkn2a-G12V"
reset_graph_directory(GRAPHS_DIR)

jun_pw <- c("KRAS", "TAB1", "MAP2K4", "MAP2K7", "MAPK8", "MAPK9", "MAPK10",
            "NR2C2", "MEN1", "JUN", "FOS", "ELK1", "JUND", "ATF2", "JDP2",
            "UBE2N", "USP28", "ZNF304", "CDT1", "KAT7",
            "CDKN2A", "CDKN2B", "TP53")


#### ---- Data preparation ---- ####

paad_cell_lines <- model_data %>%
    filter(cancer == "PAAD") %>%
    select(dep_map_id, cancer, allele) %>%
    unique()

jun_pw[!(jun_pw %in% model1_tib[model1_tib$cancer == "PAAD", ]$hugo_symbol)]


junpw_depmap <- gene_effect %>%
    inner_join(paad_cell_lines, by = c("dep_map_id")) %>%
    filter(hugo_symbol %in% jun_pw)

all(jun_pw %in% junpw_depmap$hugo_symbol)
# jun_pw[!(jun_pw %in% junpw_depmap$hugo_symbol)]


junpw_cn <- ccle_copy_number %>%
    inner_join(paad_cell_lines, by = c("dep_map_id")) %>%
    filter(hugo_symbol %in% jun_pw)

all(jun_pw %in% junpw_cn$hugo_symbol)


junpw_mut <- ccle_mutations %>%
    inner_join(paad_cell_lines, by = c("dep_map_id")) %>%
    filter(hugo_symbol %in% jun_pw) %>%
    filter(hugo_symbol != "KRAS") %>%
    select(dep_map_id, hugo_symbol, variant_classification, variant_type,
           protein_change, is_deleterious, is_cosmic_hotspot, cancer, allele) %>%
    group_by(cancer, allele, dep_map_id, hugo_symbol) %>%
    summarise(
        n_muts = n_distinct(protein_change),
        protein_change = paste(protein_change, collapse = ", "),
        is_deleterious = paste(is_deleterious, collapse = ", "),
        is_cosmic_hotspot = paste(is_cosmic_hotspot, collapse = ", ")
    ) %>%
    ungroup()

all(jun_pw %in% junpw_mut$hugo_symbol)


junpw_expr <- ccle_expression %>%
    inner_join(paad_cell_lines, by = c("dep_map_id")) %>%
    filter(hugo_symbol %in% jun_pw)

all(jun_pw %in% junpw_expr$hugo_symbol)

cn_levels <- c("het_del", "norm", "amp")

junpw_depmap %<>%
    select(cancer, dep_map_id, allele, hugo_symbol, gene_effect) %>%
    left_join(junpw_cn, by = c("hugo_symbol", "dep_map_id", "cancer", "allele")) %>%
    left_join(junpw_mut, by = c("hugo_symbol", "dep_map_id", "cancer", "allele")) %>%
    left_join(junpw_expr, by = c("hugo_symbol", "dep_map_id", "cancer", "allele")) %>%
    mutate(
        is_mutated = ifelse(!is.na(protein_change), "mut.", "not mut."),
        copy_number_label = factor(copy_number_label, levels = cn_levels)
    )

junpw_depmap %>%
    group_by(dep_map_id, hugo_symbol) %>%
    count() %>%
    filter(n > 1)


#### ---- Basic analysis ---- ####

set.seed(0)

# Box-plots of gene effect by KRAS allele
geneeffect_boxplots <- junpw_depmap %>%
    mutate(hugo_symbol = factor(hugo_symbol, levels = jun_pw)) %>%
    ggplot(aes(x = allele, y = gene_effect)) +
    facet_wrap(. ~ hugo_symbol, scales = "free", ncol = 4) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey25", size = 0.6) +
    geom_boxplot(aes(color = allele, fill = allele),
                 alpha = 0.5, outlier.shape = NA) +
    geom_jitter(aes(color = allele, shape = is_mutated),
                width = 0.25, size = 0.8) +
    scale_color_manual(values = short_allele_pal) +
    scale_fill_manual(values = short_allele_pal) +
    scale_shape_manual(values = c(17, 16)) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(strip.background = element_blank())
ggsave_wrapper(
    geneeffect_boxplots,
    plot_path(GRAPHS_DIR, "geneeffect_boxplots.svg"),
    "large"
)


# Scatter plots of gene effect vs. RNA expr for each gene
geneeffect_rnaexpr_scatter <- junpw_depmap %>%
    mutate(hugo_symbol = factor(hugo_symbol, levels = jun_pw)) %>%
    ggplot(aes(x = rna_expression, y = gene_effect)) +
    facet_wrap(. ~ hugo_symbol, scales = "free", ncol = 4) +
    geom_hline(yintercept = 0, linetype = 1, color = "grey50", size = 0.3) +
    geom_vline(xintercept = 0, linetype = 1, color = "grey50", size = 0.3) +
    geom_point(aes(color = allele,
                   shape = is_mutated,
                   size = copy_number_label)) +
    scale_color_manual(values = short_allele_pal) +
    scale_shape_manual(values = c(17, 16)) +
    scale_size_manual(values = c(1, 2, 3)) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(strip.background = element_blank())
ggsave_wrapper(
    geneeffect_rnaexpr_scatter,
    plot_path(GRAPHS_DIR, "geneeffect_rnaexpr_scatter.svg"),
    "large"
)


# JUN expression and CDKN2A gene effect.
JUNexpr_CDKN2Adep_scatter <- junpw_depmap %>%
    filter(hugo_symbol %in% c("JUN", "CDKN2A")) %>%
    select(dep_map_id, allele, hugo_symbol, gene_effect, rna_expression) %>%
    pivot_wider(c(dep_map_id, allele),
                names_from = hugo_symbol,
                values_from = c(rna_expression, gene_effect)) %>%
    ggplot(aes(x = gene_effect_JUN, rna_expression_CDKN2A)) +
    geom_point(aes(color = allele)) +
    scale_color_manual(values = short_allele_pal) +
    scale_shape_manual(values = c(17, 16)) +
    scale_size_manual(values = c(1, 2, 3)) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(strip.background = element_blank())
ggsave_wrapper(
    JUNexpr_CDKN2Adep_scatter,
    plot_path(GRAPHS_DIR, "JUNexpr_CDKN2Adep_scatter.svg"),
    "small"
)


# JUN gene effect vs. TP53 gene effect.
JUN_TP53_scatter <- junpw_depmap %>%
    filter(hugo_symbol %in% c("JUN", "TP53")) %>%
    select(dep_map_id, allele, hugo_symbol, gene_effect, rna_expression) %>%
    pivot_wider(c(dep_map_id, allele),
                names_from = hugo_symbol,
                values_from = c(rna_expression, gene_effect)) %>%
    ggplot(aes(x = gene_effect_JUN, gene_effect_TP53)) +
    geom_point(aes(color = allele)) +
    scale_color_manual(values = short_allele_pal) +
    scale_shape_manual(values = c(17, 16)) +
    scale_size_manual(values = c(1, 2, 3)) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(strip.background = element_blank())
ggsave_wrapper(
    JUN_TP53_scatter,
    plot_path(GRAPHS_DIR, "JUN_TP53_scatter.svg"),
    "small"
)

# JUN gene effect vs. MEN1 gene effect.
JUN_MEN1_scatter <- junpw_depmap %>%
    filter(hugo_symbol %in% c("JUN", "MEN1")) %>%
    select(dep_map_id, allele, hugo_symbol, gene_effect, rna_expression) %>%
    pivot_wider(c(dep_map_id, allele),
                names_from = hugo_symbol,
                values_from = c(rna_expression, gene_effect)) %>%
    ggplot(aes(x = gene_effect_JUN, gene_effect_MEN1)) +
    geom_point(aes(color = allele)) +
    scale_color_manual(values = short_allele_pal) +
    scale_shape_manual(values = c(17, 16)) +
    scale_size_manual(values = c(1, 2, 3)) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(strip.background = element_blank())
ggsave_wrapper(
    JUN_MEN1_scatter,
    plot_path(GRAPHS_DIR, "JUN_MEN1_scatter.svg"),
    "small"
)


#### ---- Comutation ---- ####

genetic_interaction_df %>%
    filter(cancer == "PAAD") %>%
    filter(hugo_symbol %in% jun_pw) %>%
    select(hugo_symbol, allele, p_val, genetic_interaction)


#### ---- JUN-regulated genes ---- ####

jun_bs <- encode_tf_bindingsites %>%
    filter(gene_set == "JUN") %>%
    unique()

jun_bs

paad_deps <- model1_tib %>%
    filter(cancer == "PAAD") %>%
    filter(!is.na(allele_aov)) %>%
    mutate(anova_pval = map_dbl(allele_aov,  ~ tidy(.x)$p.value[[1]])) %>%
    filter(anova_pval < 0.01)

jun_bs_deps <- paad_deps %>% filter(hugo_symbol %in% jun_bs$gene)

jun_bs_dep_boxplots <- jun_bs_deps %>%
    select(hugo_symbol, data) %>%
    unnest(data) %>%
    ggplot(aes(x = allele, y = gene_effect, color = allele, fill = allele)) +
    facet_wrap(~ hugo_symbol, scales = "free") +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    scale_color_manual(values = short_allele_pal) +
    scale_fill_manual(values = short_allele_pal) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(strip.background = element_blank())
ggsave_wrapper(
    jun_bs_dep_boxplots,
    plot_path(GRAPHS_DIR, "jun_bs_dep_boxplots.svg"),
    "large"
)



# Test for enrichment of JUN-regulated genes in the allele-specific
# genetic dependencies in PAAD.

all_paad_genes <- unique(model1_tib[model1_tib$cancer == "PAAD", ]$hugo_symbol)
sig_paad_genes <- unique(paad_deps$hugo_symbol)
jun_targets <- jun_bs$gene


fisher_tf_in_paad <- function(tf_genes, alternative = "two.sided") {
    tf_genes <- unlist(unique(tf_genes))
    tf_genes <- tf_genes[tf_genes %in% all_paad_genes]

    if (length(tf_genes) < 5) return(NULL)

    A <- length(setdiff(all_paad_genes, c(tf_genes, sig_paad_genes)))
    B <- length(setdiff(sig_paad_genes, tf_genes))
    C <- length(setdiff(tf_genes, sig_paad_genes))
    D <- length(intersect(sig_paad_genes, tf_genes))

    stopifnot(sum(c(A, B, C, D)) == n_distinct(all_paad_genes))

    m <- matrix(c(A, B, C, D), nrow = 2, byrow = FALSE)
    fish_res <- fisher.test(m, alternative = alternative)
    return(fish_res)
}


tf_stats <- encode_tf_bindingsites %>%
    group_by(gene_set) %>%
    nest() %>%
    ungroup() %>%
    mutate(
        fish_res = purrr::map(data, fisher_tf_in_paad),
        fish_tidy = purrr::map(fish_res, ~ tidy(.x))
    ) %>%
    unnest(fish_tidy) %>%
    janitor::clean_names()

tf_stats %>%
    mutate(p_value_adj = p.adjust(p_value, method = "BH")) %>%
    filter(p_value < 0.05 & p_value_adj < 0.25)

tf_stats_volcano <- tf_stats %>%
    mutate(
        tf_gs_len = map_dbl(data, ~ length(unlist(.x))),
        norm_log_or = log2(estimate) / sqrt(tf_gs_len),
        jun = ifelse(gene_set == "JUN", "JUN", "other"),
        labels = ifelse(p_value < 0.05 | gene_set == "JUN", gene_set, NA)
    ) %>%
    ggplot(aes(x = log2(estimate), y = -log10(p_value))) +
    geom_hline(yintercept = 0, size = 0.5, color = "grey25") +
    geom_vline(xintercept = 0, size = 0.5, color = "grey25") +
    geom_hline(yintercept = -log10(0.05),
               size = 0.7, linetype = 2, color = "tomato") +
    geom_point(aes(color = jun, size = jun), alpha = 0.6) +
    ggrepel::geom_text_repel(aes(label = labels),
                             size = 2, color = "black", family = "Arial") +
    scale_color_manual(values = c("dodgerblue", "grey50"), guide = FALSE) +
    scale_size_manual(values = c(0.9, 0.7), guide = FALSE) +
    theme_bw(base_size = 7, base_family = "Arial") +
    labs(x = "log2 OR", y = "-log10 p-value")
ggsave_wrapper(
    tf_stats_volcano,
    plot_path(GRAPHS_DIR, "tf_stats_volcano.svg"),
    "medium"
)


#### ---- JUN vs MAPK8 gene effect ---- ####

# Save scatter plot for Fig.
fignum <- 14
supp <- TRUE

# JUN vs. MAPK8 gene effect
JUN_MAPK8_scatter <- junpw_depmap %>%
    filter(hugo_symbol %in% c("JUN", "MAPK8")) %>%
    select(dep_map_id, allele, hugo_symbol, gene_effect, rna_expression) %>%
    pivot_wider(c(dep_map_id, allele),
                names_from = hugo_symbol,
                values_from = c(rna_expression, gene_effect)) %>%
    ggplot(aes(x = gene_effect_JUN, gene_effect_MAPK8)) +
    geom_point(aes(color = allele)) +
    scale_color_manual(values = short_allele_pal) +
    scale_shape_manual(values = c(17, 16)) +
    scale_size_manual(values = c(1, 2, 3)) +
    theme_bw(base_size = 7, base_family = "Arial") +
    theme(strip.background = element_blank())
ggsave_wrapper(
    JUN_MAPK8_scatter,
    plot_path(GRAPHS_DIR, "JUN_MAPK8_scatter.svg"),
    "small"
)
saveRDS(
    JUN_MAPK8_scatter,
    get_fig_proto_path("JUN_MAPK8_scatter", fignum, supp = supp)
)
