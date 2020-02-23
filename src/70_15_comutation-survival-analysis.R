# Survival analysis including genes found to comutate with KRAS alleles.

library(survminer)
library(survival)


GRAPHS_DIR <- "70_15_comutation-survival-analysis"
reset_graph_directory(GRAPHS_DIR)
reset_table_directory(GRAPHS_DIR)


# Add a column with just the KRAS allele.
cancer_coding_muts_df %<>% mutate(allele = str_remove_all(ras_allele, "KRAS_"))



# Get the survival data for given tumor samples.
get_survival_data <- function(cancers, tumor_sample_barcodes) {
    cancer_survival_df %>%
        filter(
            tumor_sample_barcode %in% cancer_coding_muts_df$tumor_sample_barcode
        ) %>%
        filter(cancer %in% !!cancers &
               tumor_sample_barcode %in% !!tumor_sample_barcodes)
}


# Get the mutation data of a cancer for specific hugo symbols.
# Optionally filter on the KRAS allele (omit "KRAS_").
get_mutation_data <- function(cancers, hugo_symbols, kras_alleles = NULL) {
    df <- cancer_coding_muts_df %>%
        filter(cancer %in% !!cancers & hugo_symbol %in% !!hugo_symbols)

    if (!is.null(kras_alleles)) df %<>% filter(allele %in% !!kras_alleles)

    return(df)
}


# Get the tumor sample barcodes for samples in a cancer.
get_cancer_sample_barcodes <- function(cancers) {
    cancer_coding_muts_df %>%
        filter(cancer %in% !!cancers) %>%
        pull(tumor_sample_barcode) %>%
        unique()
}


# Get the tumor sample barcodes for samples with the KRAS alleles in a cancer.
get_kras_allele_samples_barcodes <- function(cancers, kras_alleles) {
    cancer_coding_muts_df %>%
        filter(cancer %in% !!cancers & allele %in% !!kras_alleles) %>%
        pull(tumor_sample_barcode) %>%
        unique()
}


# A wrapper around `ggsurvplot()`.
ggsurvplot_wrapper <- function(fit, data, save_name_template,
                               tbl_col = NULL, facet_formula = NULL,
                               pal = NULL, plt_title = NULL,
                               save_size = "wide", conf_int = TRUE) {
    p <- ggsurvplot(
        fit = fit,
        data = data,
        pval = TRUE,
        conf.int = conf_int,
        risk.table = TRUE,
        risk.table.col = tbl_col,
        surv.median.line = "hv",
        fontsize = 3,
        font.family = "arial",
        palette = pal
    ) +
        ggtitle(plt_title)
    p <- style_ggsurvminer_plot(p)

    if (!is.null(facet_formula)) {
        p$plot <- p$plot + facet_grid(facet_formula)
    }

    surv_plot <- patch_ggsurvplot(p)
    ggsave_wrapper(
        surv_plot,
        plot_path(GRAPHS_DIR, glue(save_name_template)),
        save_size
    )
}


# Model the effect of comutation of another gene in the samples with a
# specific KRAS allele.
alleleonly_samples_comutation_sa <- function(cancer,
                                             allele,
                                             hugo_symbol,
                                             genetic_interaction,
                                             min_samples_per_group = 5,
                                             p_val = 0.05) {
    mut_data <- get_mutation_data(cancer, hugo_symbol, allele)
    comut_samples <- mut_data$tumor_sample_barcode
    all_allele_samples <- get_kras_allele_samples_barcodes(cancer, allele)
    surv_data <- get_survival_data(cancer, all_allele_samples)

    surv_data %<>%
        mutate(comutation = tumor_sample_barcode %in% comut_samples)

    # Return early if not enough samples to be worth it.
    if (sum(surv_data$comutation) < min_samples_per_group) return()
    if (sum(!surv_data$comutation) < min_samples_per_group) return()

    fit <- survfit(Surv(time = time, event = status) ~ comutation,
                   data = surv_data)
    fit_diff <- survdiff(Surv(time = time, event = status) ~ comutation,
                         data = surv_data)

    title <- glue("{cancer} - KRAS {allele} - comutation with {hugo_symbol}")
    fname <- glue("samplesonly_{cancer}_{allele}.txt")
    write_survfit_summary(fit, title = title, fname = fname)
    write_survdiff_summary(fit_diff, title = "", fname = fname)

    # Return early if p-value of `survdiff` is above cut-off
    if (survdiff_pval(fit_diff) >= p_val) return()

    plt_title <- glue(
        "{cancer} - KRAS {allele} (only)
        comutation with {hugo_symbol} ({genetic_interaction})"
    )
    ggsurvplot_wrapper(
        fit, surv_data,
        glue("samplesonly_{cancer}-{allele}-{hugo_symbol}.svg"),
        tbl_col = "comutation",
        pal = c("grey50", "mediumpurple2"),
        plt_title = plt_title
    )

}


# Model survival analysis on the KRAS allele and comutation of another gene
# only in samples with a mutant KRAS.
krasallele_comutation_sa <- function(cancer,
                                     allele,
                                     hugo_symbol,
                                     genetic_interaction,
                                     min_samples_per_group = 5,
                                     p_val = 0.05) {
    mut_data <- get_mutation_data(cancer, hugo_symbol)
    comut_samples <- mut_data$tumor_sample_barcode
    cancer_samples <- get_cancer_sample_barcodes(cancer)
    surv_data <- get_survival_data(cancer, cancer_samples)

    surv_data %<>% mutate(
        comutation = tumor_sample_barcode %in% comut_samples,
        krasallele = kras_allele == !!allele
    )

    # Remove groups with fewer than `min_samples_per_group`.
    surv_data %<>%
        group_by(comutation, krasallele) %>%
        filter(n_distinct(tumor_sample_barcode) >= !!min_samples_per_group) %>%
        ungroup()

    # If there is no data, return early.
    if (nrow(surv_data) == 0) return()

    # If fewer than 3 groups, return early
    n_groups <- surv_data %>%
        group_by(comutation, krasallele) %>%
        count() %>%
        nrow()
    if (n_groups < 3) return()

    fit <- survfit(
        Surv(time = time, event = status) ~ krasallele + comutation,
        data = surv_data
    )

    fit_coxph <- coxph(
        Surv(time = time, event = status) ~ krasallele + comutation,
        data = surv_data
    )

    title <- glue("{cancer} - KRAS {allele} - comutation with {hugo_symbol}")
    fname <- glue("allsamples_{cancer}_{allele}.txt")
    write_survfit_summary(fit, title = title, fname = fname)
    write_coxph_summary(fit_coxph, title = "", fname = fname)

    # Return early if p-value of `coxph` is above cut-off
    if (coxph_logtest_pval(fit_coxph) >= p_val) return()

    plt_title <- glue(
        "{cancer} - KRAS {allele}
        comutation with {hugo_symbol} ({genetic_interaction})"
    )

    # Weird way of making a palette with the correct colors
    pal <- tibble::tribble(
        ~ comutation, ~krasallele, ~color,
        FALSE, FALSE, "grey50",
        FALSE, TRUE, "grey25",
        TRUE, FALSE, "plum2",
        TRUE, TRUE, "mediumpurple2"
    ) %>%
        inner_join(surv_data, by = c("comutation", "krasallele")) %>%
        select(comutation, krasallele, color) %>%
        unique() %>%
        pull(color)

    ggsurvplot_wrapper(
        fit, surv_data,
        glue("allsamples_{cancer}-{allele}-{hugo_symbol}.svg"),
        pal = pal,
        plt_title = plt_title,
        conf_int = FALSE
    )
}


# Model survival analysis on the KRAS allele and comutation of another gene
# only using samples with a KRAS mutation.
krasmutsamples_krasallele_comutation_sa <- function(cancer,
                                                    allele,
                                                    hugo_symbol,
                                                    genetic_interaction,
                                                    min_samples_per_group = 5,
                                                    p_val = 0.05) {
    mut_data <- get_mutation_data(cancer, hugo_symbol)
    comut_samples <- mut_data$tumor_sample_barcode
    surv_data <- get_survival_data(cancer, get_cancer_sample_barcodes(cancer))

    surv_data %<>%
        filter(kras_allele != "WT") %>%
        mutate(
            comutation = tumor_sample_barcode %in% comut_samples,
            krasallele = kras_allele == !!allele
        )

    # Remove groups with fewer than `min_samples_per_group`.
    surv_data %<>%
        group_by(comutation, krasallele) %>%
        filter(n_distinct(tumor_sample_barcode) >= !!min_samples_per_group) %>%
        ungroup()

    # If there is no data, return early.
    if (nrow(surv_data) == 0) return()

    # If fewer than 3 groups, return early
    n_groups <- surv_data %>%
        group_by(comutation, krasallele) %>%
        count() %>%
        nrow()
    if (n_groups < 3) return()

    fit <- survfit(
        Surv(time = time, event = status) ~ krasallele + comutation,
        data = surv_data
    )

    fit_coxph <- coxph(
        Surv(time = time, event = status) ~ krasallele + comutation,
        data = surv_data
    )

    title <- glue("{cancer} - KRAS {allele} - comutation with {hugo_symbol}")
    fname <- glue("krasmutsamples_{cancer}_{allele}.txt")
    write_survfit_summary(fit, title = title, fname = fname)
    write_coxph_summary(fit_coxph, title = "", fname = fname)

    # Return early if p-value of `coxph` is above cut-off
    if (coxph_logtest_pval(fit_coxph) >= p_val) return()

    plt_title <- glue(
        "{cancer} - KRAS {allele} (KRAS muts. only)
        comutation with {hugo_symbol} ({genetic_interaction})"
    )

    # Weird way of making a palette with the correct colors
    pal <- tibble::tribble(
        ~ comutation, ~krasallele, ~color,
        FALSE, FALSE, "grey50",
        FALSE, TRUE, "grey25",
        TRUE, FALSE, "plum2",
        TRUE, TRUE, "mediumpurple2"
    ) %>%
        inner_join(surv_data, by = c("comutation", "krasallele")) %>%
        select(comutation, krasallele, color) %>%
        unique() %>%
        pull(color)

    ggsurvplot_wrapper(
        fit, surv_data,
        glue("krasmutsamples_{cancer}-{allele}-{hugo_symbol}.svg"),
        pal = pal,
        plt_title = plt_title,
        conf_int = FALSE
    )
}


# Model survival analysis on the KRAS allele and comutation of another gene
# only using samples with the KRAS allele or WT KRAS.
alleleorwt_krasallele_comutation_sa <- function(cancer,
                                                    allele,
                                                    hugo_symbol,
                                                    genetic_interaction,
                                                    min_samples_per_group = 5,
                                                    p_val = 0.05) {
    mut_data <- get_mutation_data(cancer, hugo_symbol)
    comut_samples <- mut_data$tumor_sample_barcode
    surv_data <- get_survival_data(cancer, get_cancer_sample_barcodes(cancer))

    surv_data %<>%
        filter(kras_allele %in% c("WT", !!allele)) %>%
        mutate(
            comutation = tumor_sample_barcode %in% comut_samples,
            krasallele = kras_allele == !!allele
        )

    # Remove groups with fewer than `min_samples_per_group`.
    surv_data %<>%
        group_by(comutation, krasallele) %>%
        filter(n_distinct(tumor_sample_barcode) >= !!min_samples_per_group) %>%
        ungroup()

    # If there is no data, return early.
    if (nrow(surv_data) == 0) return()

    # If fewer than 3 groups, return early
    n_groups <- surv_data %>%
        group_by(comutation, krasallele) %>%
        count() %>%
        nrow()
    if (n_groups < 3) return()

    fit <- survfit(
        Surv(time = time, event = status) ~ krasallele + comutation,
        data = surv_data
    )

    fit_coxph <- coxph(
        Surv(time = time, event = status) ~ krasallele + comutation,
        data = surv_data
    )

    title <- glue("{cancer} - KRAS {allele} - comutation with {hugo_symbol}")
    fname <- glue("alleleorwt_{cancer}_{allele}.txt")
    write_survfit_summary(fit, title = title, fname = fname)
    write_coxph_summary(fit_coxph, title = "", fname = fname)

    # Return early if p-value of `coxph` is above cut-off
    if (coxph_logtest_pval(fit_coxph) >= p_val) return()

    plt_title <- glue(
        "{cancer} - KRAS {allele} (KRAS allele or KRAS WT only)
        comutation with {hugo_symbol} ({genetic_interaction})"
    )

    # Weird way of making a palette with the correct colors
    pal <- tibble::tribble(
        ~ comutation, ~krasallele, ~color,
        FALSE, FALSE, "grey50",
        FALSE, TRUE, "grey25",
        TRUE, FALSE, "plum2",
        TRUE, TRUE, "mediumpurple2"
    ) %>%
        inner_join(surv_data, by = c("comutation", "krasallele")) %>%
        select(comutation, krasallele, color) %>%
        unique() %>%
        pull(color)

    ggsurvplot_wrapper(
        fit, surv_data,
        glue("alleleorwt_{cancer}-{allele}-{hugo_symbol}.svg"),
        pal = pal,
        plt_title = plt_title,
        conf_int = FALSE
    )
}



genetic_interaction_df %>%
    select(cancer, allele, hugo_symbol, genetic_interaction) %>%
    unique() %>%
    arrange(cancer, allele, hugo_symbol) %T>%
    pwalk(alleleonly_samples_comutation_sa) %T>%
    pwalk(krasallele_comutation_sa) %T>%
    pwalk(krasmutsamples_krasallele_comutation_sa) %T>%
    pwalk(alleleorwt_krasallele_comutation_sa)
