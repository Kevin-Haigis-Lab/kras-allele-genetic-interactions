# Comutation analysis of the alleles of other oncogenes with KRAS alleles.

GRAPHS_DIR <- "20_75_allele-to-allele-comutation-analysis"
reset_graph_directory(GRAPHS_DIR)


oncogenes_to_test <- cosmic_cgc_df %>%
    distinct(cancer, hugo_symbol, name, tier, role_in_cancer, mutation_types) %>%
    filter(cancer %in% c("COAD", "LUAD", "PAAD", "MM")) %>%
    filter(hugo_symbol != "KRAS")
oncogenes_to_test


panel_gene_lists <- cancer_full_coding_muts_df %>%
    distinct(cancer, hugo_symbol, dataset, target)

datasets_for_gene <- function(cancer, gene) {
    panel_gene_lists %>%
        filter(hugo_symbol == !!gene & cancer == !!cancer) %>%
        u_pull(dataset)
}
datasets_for_gene <- memoise::memoise(datasets_for_gene)


kras_alleles_to_test_against <- genetic_interaction_df %>%
    distinct(cancer, kras_allele)

kras_alleles_for_cancer <- function(cancer) {
    kras_alleles_to_test_against %>%
        filter(cancer == !!cancer) %>%
        pull(kras_allele)
}

make_lgl_table <- function(x, y) {
    f <- function(a) { factor(a, levels = c("FALSE", "TRUE")) }
    table(f(x), f(y))
}

tidy_ct_table <- function(ct_tbl) {
    as.data.frame(ct_tbl) %>%
        as_tibble() %>%
        rename(kras = Var1, gene = Var2)
}


test_for_allele_comutation_per_kras_allele <- function(comut_df, cancer) {
    test_comutation <- function(ras_allele) {
        ct_tbl <- comut_df %>%
            mutate(is_allele = ras_allele == !!ras_allele) %$%
            make_lgl_table(is_allele, is_allele_mutant)

        tidy_tbl <- tidy_ct_table(ct_tbl)

        fisher.test(ct_tbl) %>%
            broom::tidy() %>%
            bind_cols(tibble(cont_tbl = list(tidy_tbl)))
    }

    tibble(ras_allele = kras_alleles_for_cancer(cancer)) %>%
        mutate(comut_res = map(ras_allele, test_comutation)) %>%
        unnest(comut_res) %>%
        janitor::clean_names()

}


allele_to_allele_comutation_analysis <- function(cancer,
                                                 hugo_symbol,
                                                 amino_acid_change,
                                                 ...) {
    datasets_with_gene <- datasets_for_gene(cancer, hugo_symbol)
    comutation_df <- cancer_full_coding_muts_df %>%
        filter(dataset %in% !!datasets_with_gene) %>%
        group_by(tumor_sample_barcode) %>%
        summarise(
            ras_allele = unique(ras_allele),
            is_allele_mutant = any(hugo_symbol == !!hugo_symbol &
                                   amino_acid_change == !!amino_acid_change)
        ) %>%
        test_for_allele_comutation_per_kras_allele(cancer = cancer) %>%
        mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
        select(-method, -alternative)

}


allele_to_allele_comutation_res <- cancer_full_coding_muts_df %>%
    filter(cancer != "SKCM") %>%
    inner_join(oncogenes_to_test %>% distinct(cancer, hugo_symbol, role_in_cancer),
               by = c("cancer", "hugo_symbol")) %>%
    mutate(amino_acid_change = ifelse(amino_acid_change == "",
                                      mutation_type,
                                      amino_acid_change)) %>%
    count(cancer, hugo_symbol, role_in_cancer, amino_acid_change,
          name = "allele_ct") %>%
    filter(allele_ct > 10)

nrow(allele_to_allele_comutation_res)


allele_to_allele_comutation_res$comut_res <- pmap(
    allele_to_allele_comutation_res,
    allele_to_allele_comutation_analysis
)


allele_to_allele_comutation_res %>%
    unnest(comut_res) %>%
    filter(adj_p_value < 0.2) %>%
    select(-cont_tbl) %>%
    filter(str_detect(role_in_cancer, "oncogene")) %>%
    arrange(cancer, hugo_symbol, allele_ct) %>%
    knitr::kable(format = "markdown", digits = 3)

#> |cancer |hugo_symbol |role_in_cancer        |amino_acid_change | allele_ct|ras_allele | estimate| p_value| conf_low| conf_high| adj_p_value|
#> |:------|:-----------|:---------------------|:-----------------|---------:|:----------|--------:|-------:|--------:|---------:|-----------:|
#> |COAD   |BRAF        |oncogene, fusion      |V600E             |       473|KRAS_G12C  |    0.000|   0.000|    0.000|     0.279|       0.000|
#> |COAD   |BRAF        |oncogene, fusion      |V600E             |       473|KRAS_G12D  |    0.028|   0.000|    0.003|     0.102|       0.000|
#> |COAD   |BRAF        |oncogene, fusion      |V600E             |       473|KRAS_G12S  |    0.000|   0.000|    0.000|     0.426|       0.001|
#> |COAD   |BRAF        |oncogene, fusion      |V600E             |       473|KRAS_G12V  |    0.000|   0.000|    0.000|     0.077|       0.000|
#> |COAD   |BRAF        |oncogene, fusion      |V600E             |       473|KRAS_G13D  |    0.000|   0.000|    0.000|     0.089|       0.000|
#> |COAD   |BRAF        |oncogene, fusion      |V600E             |       473|KRAS_A146T |    0.000|   0.000|    0.000|     0.263|       0.000|
#> |COAD   |BRAF        |oncogene, fusion      |V600E             |       473|KRAS_G12A  |    0.000|   0.000|    0.000|     0.378|       0.000|
#> |COAD   |BRAF        |oncogene, fusion      |V600E             |       473|KRAS_Q61H  |    0.000|   0.069|    0.000|     1.106|       0.087|
#> |COAD   |BRAF        |oncogene, fusion      |V600E             |       473|KRAS_A146V |    0.000|   0.104|    0.000|     1.371|       0.116|
#> |COAD   |CTNNB1      |oncogene, fusion      |S45F              |        39|KRAS_A146T |    4.423|   0.017|    1.123|    12.683|       0.174|
#> |COAD   |PIK3CA      |oncogene              |E545G             |        21|KRAS_G12D  |    3.655|   0.009|    1.243|     9.726|       0.088|
#> |COAD   |PIK3CA      |oncogene              |E545G             |        21|KRAS_A146T |    6.166|   0.018|    1.149|    21.515|       0.088|
#> |COAD   |PIK3CA      |oncogene              |C420R             |        37|KRAS_G13D  |    4.250|   0.001|    1.744|     9.403|       0.009|
#> |COAD   |PIK3CA      |oncogene              |Q546K             |        42|KRAS_G12D  |    3.688|   0.000|    1.783|     7.293|       0.003|
#> |COAD   |PIK3CA      |oncogene              |Q546K             |        42|KRAS_G12V  |    2.579|   0.022|    1.024|     5.724|       0.111|
#> |COAD   |PIK3CA      |oncogene              |E545K             |       244|KRAS_G12D  |    2.079|   0.000|    1.480|     2.881|       0.000|
#> |COAD   |PIK3CA      |oncogene              |E545K             |       244|KRAS_G12V  |    2.049|   0.000|    1.381|     2.972|       0.002|
#> |COAD   |PIK3CA      |oncogene              |E545K             |       244|KRAS_G13D  |    2.074|   0.001|    1.367|     3.064|       0.002|
#> |COAD   |TCF7L2      |oncogene, fusion      |R488C             |        28|KRAS_G12V  |    3.339|   0.017|    1.094|     8.637|       0.173|
#> |COAD   |TP53        |oncogene, TSG, fusion |C176Y             |        16|KRAS_G12V  |    4.941|   0.009|    1.339|    15.516|       0.086|
#> |COAD   |TP53        |oncogene, TSG, fusion |X307_splice       |        21|KRAS_A146T |    6.133|   0.018|    1.143|    21.400|       0.178|
#> |COAD   |TP53        |oncogene, TSG, fusion |R342*             |        57|KRAS_Q61L  |   10.534|   0.019|    1.152|    46.169|       0.193|
#> |LUAD   |BRAF        |oncogene, fusion      |V600E             |        76|KRAS_G12C  |    0.160|   0.001|    0.019|     0.600|       0.012|
#> |LUAD   |BRAF        |oncogene, fusion      |V600E             |        76|KRAS_G12D  |    0.000|   0.033|    0.000|     0.910|       0.125|
#> |LUAD   |BRAF        |oncogene, fusion      |V600E             |        76|KRAS_G12V  |    0.174|   0.042|    0.004|     1.007|       0.125|
#> |LUAD   |EGFR        |oncogene              |L747_P753delinsS  |        32|KRAS_G12C  |    0.000|   0.010|    0.000|     0.724|       0.094|
#> |LUAD   |EGFR        |oncogene              |S768I             |        44|KRAS_G12C  |    0.000|   0.003|    0.000|     0.557|       0.025|
#> |LUAD   |EGFR        |oncogene              |T790M             |       120|KRAS_G12A  |    0.000|   0.051|    0.000|     1.012|       0.115|
#> |LUAD   |EGFR        |oncogene              |T790M             |       120|KRAS_G12C  |    0.100|   0.000|    0.012|     0.372|       0.000|
#> |LUAD   |EGFR        |oncogene              |T790M             |       120|KRAS_G12D  |    0.000|   0.003|    0.000|     0.569|       0.008|
#> |LUAD   |EGFR        |oncogene              |T790M             |       120|KRAS_G12V  |    0.000|   0.000|    0.000|     0.408|       0.002|
#> |LUAD   |EGFR        |oncogene              |E746_A750del      |       315|KRAS_G12A  |    0.000|   0.000|    0.000|     0.383|       0.000|
#> |LUAD   |EGFR        |oncogene              |E746_A750del      |       315|KRAS_G12C  |    0.000|   0.000|    0.000|     0.069|       0.000|
#> |LUAD   |EGFR        |oncogene              |E746_A750del      |       315|KRAS_G12D  |    0.058|   0.000|    0.001|     0.329|       0.000|
#> |LUAD   |EGFR        |oncogene              |E746_A750del      |       315|KRAS_G12V  |    0.000|   0.000|    0.000|     0.154|       0.000|
#> |LUAD   |EGFR        |oncogene              |E746_A750del      |       315|KRAS_G13C  |    0.000|   0.033|    0.000|     0.893|       0.059|
#> |LUAD   |EGFR        |oncogene              |E746_A750del      |       315|KRAS_G13D  |    0.000|   0.111|    0.000|     1.354|       0.166|
#> |LUAD   |EGFR        |oncogene              |L858R             |       356|KRAS_G12A  |    0.000|   0.000|    0.000|     0.329|       0.000|
#> |LUAD   |EGFR        |oncogene              |L858R             |       356|KRAS_G12C  |    0.049|   0.000|    0.010|     0.144|       0.000|
#> |LUAD   |EGFR        |oncogene              |L858R             |       356|KRAS_G12D  |    0.100|   0.000|    0.012|     0.369|       0.000|
#> |LUAD   |EGFR        |oncogene              |L858R             |       356|KRAS_G12V  |    0.000|   0.000|    0.000|     0.132|       0.000|
#> |LUAD   |EGFR        |oncogene              |L858R             |       356|KRAS_G13C  |    0.000|   0.014|    0.000|     0.767|       0.025|
#> |LUAD   |EGFR        |oncogene              |L858R             |       356|KRAS_G13D  |    0.000|   0.072|    0.000|     1.164|       0.108|
#> |LUAD   |ERBB2       |oncogene, fusion      |Y772_A775dup      |        57|KRAS_G12C  |    0.000|   0.000|    0.000|     0.393|       0.002|
#> |LUAD   |ERBB2       |oncogene, fusion      |Y772_A775dup      |        57|KRAS_G12V  |    0.000|   0.032|    0.000|     0.890|       0.142|
#> |LUAD   |TP53        |oncogene, TSG, fusion |C242F             |        13|KRAS_Q61L  |   47.167|   0.001|    4.742|   241.754|       0.012|
#> |LUAD   |TP53        |oncogene, TSG, fusion |R249M             |        18|KRAS_G12V  |    5.087|   0.007|    1.412|    15.315|       0.062|


##> COAD - CTNNB1
allele_to_allele_comutation_res %>%
    filter(cancer == "COAD" & hugo_symbol == "CTNNB1") %>%
    unnest(comut_res) %>%
    select(amino_acid_change, ras_allele, p_value, estimate) %>%
    knitr::kable(format = "markdown", digits = 3)


##> COAD - TCF7L2
genetic_interaction_df %>%
    filter(cancer == "COAD" & hugo_symbol == "TCF7L2")

cancer_full_coding_muts_df %>%
    filter(cancer == "COAD" & hugo_symbol == "TCF7L2") %>%
    count(amino_acid_change, sort = TRUE)

allele_to_allele_comutation_res %>%
    filter(cancer == "COAD" & hugo_symbol == "TCF7L2") %>%
    unnest(comut_res) %>%
    select(amino_acid_change, ras_allele, p_value, estimate) %>%
    knitr::kable(format = "markdown", digits = 3)


#### ---- Volcano plots of oncogene comutation results ---- ####

plt_data <- allele_to_allele_comutation_res %>%
    unnest(comut_res) %>%
    filter(str_detect(role_in_cancer, "oncogene")) %>%
    mutate(num_samples = map_dbl(cont_tbl, ~ sum(.x$Freq)),
           min_estimate = 1/num_samples,
           estimate = ifelse(estimate > 0, estimate, min_estimate),
           log_or = log(estimate),
           log_adj_p_value = -log10(adj_p_value),
           sig = adj_p_value < 0.2,
           allele = str_remove(ras_allele, "KRAS_"),
           label = glue("{hugo_symbol} {amino_acid_change} - {allele}"),
           label = ifelse(sig, label, NA_character_))

comutation_volcano_plot <- function(df, xlims) {
    p <- df %>% ggplot(aes(x = log_or, y = log_adj_p_value))

    if (min(df$log_adj_p_value) < 1) {
        p <- p +
            geom_hline(yintercept = -log10(0.2), lty = 2, size = 0.6)
    }

    p +
        geom_vline(xintercept = 0, lty = 2, size = 0.6) +
        geom_point(aes(color = sig), alpha = 0.6, size = 2) +
        ggrepel::geom_label_repel(
            aes(label = label),
            size = 2,
            label.padding = unit(0.8, "mm"),
            point.padding = 1e-06,
            label.r = unit(0.2, "mm"),
            label.size = unit(0.1, "mm"),
            min.segment.length = 0.5,
            family = "Arial"
        ) +
        scale_x_continuous(limits = xlims,
                           expand = expansion(mult = c(0.02, 0.02))) +
        scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "tomato"),
                           labels = c("FALSE" = "adj. p-value ≥ 0.2",
                                      "TRUE" = "adj. p-value < 0.2")) +
        theme_bw(base_size = 7, base_family = "Arial") +
        labs(x = "log odds ratio",
             y = "adj. p-value (-log10)",
             color = NULL)
}

stacked_comutation_volcano_plot <- function(df, cancer, ycut) {
    xlims <- c(min(df$log_or), max(df$log_or)) + c(-1, 1)
    plt_bottom <- df %>%
        filter(log_adj_p_value < ycut) %>%
        comutation_volcano_plot(xlims) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.02)))
    plt_top <- df %>%
        filter(log_adj_p_value > ycut) %>%
        comutation_volcano_plot(xlims) +
        scale_y_continuous(expand = expansion(mult = c(0.04, 0.04))) +
        theme(axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none",
              plot.title = element_text(hjust = 0.5)) +
        labs(title = cancer)
    p <- (plt_top / plt_bottom) +
        plot_layout(heights = c(1, 2))
    return(p)
}


plot_and_save_stacked_comut_volcano <- function(cancer, ycut, data) {
    p <- data %>%
        filter(cancer == !!cancer) %>%
        stacked_comutation_volcano_plot(cancer = cancer, ycut = ycut)
    ggsave_wrapper(
        p,
        plot_path(GRAPHS_DIR, glue("{cancer}_comutation-volcano.svg")),
        "medium"
    )
}

a <- tibble(cancer = c("COAD", "LUAD", "PAAD", "MM"),
            ycut = c(10, 6, 6, 6)) %>%
    pmap(plot_and_save_stacked_comut_volcano, data = plt_data)


pik3ca_volcano <- plt_data %>%
    filter(cancer == "COAD" & hugo_symbol == "PIK3CA") %>%
    comutation_volcano_plot(xlims = c(NA, NA)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    theme(plot.title = element_markdown(hjust = 0.5)) +
    labs(title = "COAD - *PIK3CA* allele-specific comutation with *KRAS*")
ggsave_wrapper(pik3ca_volcano,
               plot_path(GRAPHS_DIR, "COAD_PIK3CA_comutation-volcano.svg"),
               "medium")



# Heatmap of PIK3CA allele-specific comutation with KRAS
y_text_glue <- "{amino_acid_change} <span style='color:#8c8c8c'>({allele_ct})</span>"

pik3ca_comut_heatmap <- allele_to_allele_comutation_res %>%
    unnest(comut_res) %>%
    filter(cancer == "COAD" & hugo_symbol == "PIK3CA") %>%
    mutate(
        log_or = ifelse(estimate == 0, NA, log(estimate)),
        log_adj_p_value = p.adjust(adj_p_value, method = "BH"),
        pik3ca_codon = str_extract(amino_acid_change, "[:digit:]+"),
        pik3ca_codon = as.numeric(pik3ca_codon),
        ras_allele = str_remove(ras_allele, "KRAS_"),
        ras_allele = factor_alleles(ras_allele),
        amino_acid_change_allele_ct = glue(y_text_glue),
        amino_acid_change_allele_ct = fct_reorder(amino_acid_change_allele_ct,
                                                  pik3ca_codon),
        stars = assign_stars(adj_p_value),
        pik3ca_codon = fct_rev(factor(pik3ca_codon))
    ) %>%
    ggplot(aes(ras_allele, amino_acid_change_allele_ct)) +
    facet_grid(pik3ca_codon ~ ., scales = "free_y", space = "free_y") +
    geom_tile(aes(fill = log_or), color = "white") +
    geom_text(aes(label = stars), size = 2) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradient2(low = "#0CA2B3", mid = "white", high = "#B34C00",
                         midpoint = 0, na.value = "grey95") +
    theme_bw(base_size = 11, base_family = "Arial") +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown(),
          axis.text.x = element_markdown(),
          axis.text.y = element_markdown(),
          axis.ticks = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.spacing = unit(0.5, "mm")) +
    labs(x = "*KRAS* alleles",
         y = "*PIK3CA* alleles",
         fill = "log odds\nof comutation")
ggsave_wrapper(pik3ca_comut_heatmap,
               plot_path(GRAPHS_DIR, "COAD_PIK3CA_comutation-heatmap.svg"),
               "medium")
