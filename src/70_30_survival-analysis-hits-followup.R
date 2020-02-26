# Follow-up and further inspection of the samples with interesting results
# from the survival analysis.

# `survival_analysis_hits`: tibble from "lib/survival_analysis_hits.R"
#       `cancer_oncogenes`: tibble from "lib/cancer-specific_helpers.R"

GRAPHS_DIR <- "70_30_survival-analysis-hits-followup"
reset_graph_directory(GRAPHS_DIR)
reset_table_directory(GRAPHS_DIR)


#### ---- Prepare tumor sample data frame ---- ####

# Data frame of samples with mutations in survival analysis hit genes.
sahits_sample_df <- cancer_coding_muts_df %>%
    mutate(allele = str_remove_all(ras_allele, "KRAS_")) %>%
    right_join(survival_analysis_hits,
               by = c("cancer", "hugo_symbol")) %>%
    select(cancer, allele, tumor_sample_barcode, hugo_symbol,
           interaction_allele, mutation_type, VAF, comutation_interaction) %>%
    unique()





get_other_oncogene_mutations <- function(tsb, cancer) {
    oncogenes <-  c(cancer_oncogenes[[cancer]], "KRAS")
    cancer_coding_muts_df %>%
        filter(tumor_sample_barcode %in% !!tsb &
               hugo_symbol %in% !!oncogenes) %>%
        pull(hugo_symbol) %>%
        unlist() %>%
        unique() %>%
        sort() %>%
        paste(collapse = ", ")
}
get_other_oncogene_mutations <- memoise::memoise(get_other_oncogene_mutations)


get_oncogene_label_order <- function(lbls) {
    lbls <- unique(lbls)
    lbls <- lbls[lbls != ""]
    new_lvls <- fct_reorder2(lbls, order(lbls), -str_count(lbls, ",")) %>%
        levels()
    return(c("none", new_lvls))
}


barplot_of_oncogene_comuts <- function(cancer, hugo_symbol,
                                       interaction_allele, data) {
    oncogene_lbl_order <- get_oncogene_label_order(data$other_oncogene_muts)

    p <- data %>%
        mutate(
            has_interaction_allele = allele == !!interaction_allele,
            other_oncogene_muts = ifelse(
                other_oncogene_muts == "", "none", other_oncogene_muts
            )
        ) %>%
        group_by(has_interaction_allele, other_oncogene_muts) %>%
        summarise(num = n_distinct(tumor_sample_barcode)) %>%
        ungroup() %>%
        tidyr::complete(
            has_interaction_allele, other_oncogene_muts, fill = list(num = 0)
        ) %>%
        mutate(other_oncogene_muts = factor(other_oncogene_muts,
                                            levels = oncogene_lbl_order)) %>%
        ggplot(aes(x = other_oncogene_muts, y = num)) +
        geom_col(aes(fill = has_interaction_allele), position = "dodge") +
        scale_fill_manual(values = c("grey60", "grey30")) +
        scale_y_continuous(expand = expand_scale(mult = c(0, 0.02))) +
        theme_bw(base_size = 7, base_family = "arial") +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 30, hjust = 1),
            axis.ticks = element_blank()
        ) +
        labs(
            x = "oncogenes mutated",
            y = "num. tumor samples",
            title = glue(
                "{cancer} - samples with a mutation in {hugo_symbol}
                separated by those with KRAS {interaction_allele}"),
            fill = glue("has {interaction_allele}")
        )

    ggsave_wrapper(
        p,
        plot_path(GRAPHS_DIR, glue("other-oncogene-muts_{hugo_symbol}_{interaction_allele}_{cancer}.svg")),
        "small"
    )
}



plot_tmb_of_oncogene_comuts <- function(cancer,
                                        interaction_allele,
                                        data) {
    tsb_allele <- data %>%
        select(tumor_sample_barcode, allele, hugo_symbol) %>%
        unique() %>%
        mutate(has_interaction_allele = allele == !!interaction_allele)

    p <- cancer_coding_muts_df %>%
        filter(tumor_sample_barcode %in% data$tumor_sample_barcode) %>%
        group_by(tumor_sample_barcode) %>%
        summarise(tmb = n()) %>%
        ungroup() %>%
        right_join(tsb_allele, by = "tumor_sample_barcode") %>%
        ggplot(aes(x = hugo_symbol, y = tmb)) +
        ggbeeswarm::geom_quasirandom(aes(color = has_interaction_allele),
                                     size = 1, alpha = 0.8) +
        scale_color_manual(values = c("grey60", "grey20")) +
        scale_y_continuous(trans = "log10") +
        theme_bw(base_size = 7, base_family = "arial") +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.ticks = element_blank(),
            axis.title.x = element_blank()
        ) +
        labs(
            y = "tumor mutational burden (coding mutations)",
            title = glue("{cancer} - TMB in samples with comutation with {interaction_allele}"),
            color = glue("has {interaction_allele}")
        )
    ggsave_wrapper(
        p,
        plot_path(GRAPHS_DIR, glue("tmb-comut-samples_{interaction_allele}_{cancer}.svg")),
        "small"
    )
}


sahits_sample_df %>%
    mutate(other_oncogene_muts = map2_chr(
        tumor_sample_barcode, cancer, get_other_oncogene_mutations
    )) %>%
    group_by(cancer, hugo_symbol, interaction_allele) %>%
    nest() %T>%
    pwalk(barplot_of_oncogene_comuts) %>%
    unnest(data) %>%
    group_by(cancer, interaction_allele) %>%
    nest() %T>%
    pwalk(plot_tmb_of_oncogene_comuts)
