# Investigating the relationship between JUN and KRAS G12V in PAAD.

GRAPHS_DIR <- "90_25_paad_depmap_jun-cdkn2a-G12V"
reset_graph_directory(GRAPHS_DIR)

jun_pw <- c("KRAS", "MAP2K4", "MAP2K7", "MAPK8", "MAPK9", "MAPK10", "MEN1",
            "JUN", "FOS", "ELK1", "JUND", "ATF2", "JDP2",
            "USP28", "ZNF304", "CDT1", "KAT7",
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
           protein_change, is_deleterious, is_cosmic_hotspot, cancer, allele)

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
