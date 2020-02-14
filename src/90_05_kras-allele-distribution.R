# Plot the distribution of KRAS mutations.

GRAPHS_DIR <- "90_05_kras-allele-distribution"
reset_graph_directory(GRAPHS_DIR)
reset_table_directory(GRAPHS_DIR)


# A data frame of the KRAS allele for each `tumor_sample_barcode`.
alleles_df <- cancer_full_coding_muts_df %>%
    group_by(tumor_sample_barcode) %>%
    slice(1) %>%
    ungroup() %>%
    select(cancer, tumor_sample_barcode, dataset, target,
           is_hypermutant, ras, ras_allele) %>%
    unique()

# The alleles to use for plotting.
alleles_to_keep <- names(short_allele_pal)
alleles_to_keep <- alleles_to_keep[alleles_to_keep != "Other"]

# A data frame of the distribution of alleles across cancer.
get_allele_distribtion_dataframe <- function(with_other = TRUE) {
    df <- alleles_df %>%
        filter(!is_hypermutant) %>%
        mutate(ras_allele = str_remove(ras_allele, "KRAS_"))

    if (with_other) {
        df %<>%
            mutate(
                ras_allele = fct_other(ras_allele, keep = alleles_to_keep)
            ) %>%
            group_by(cancer) %>%
            mutate(
                ras_allele = as.character(fct_lump(ras_allele, prop = 0.01)),
            )
    }
    df %<>%
        group_by(cancer) %>%
        mutate(
            num_cancer_samples = n_distinct(tumor_sample_barcode)
        ) %>%
        group_by(cancer, ras, ras_allele, num_cancer_samples) %>%
        summarise(num_allele_samples = n_distinct(tumor_sample_barcode)) %>%
        ungroup()

    df %<>%
        filter(cancer != "SKCM") %>%
        mutate(allele_freq = num_allele_samples / num_cancer_samples)

    return(df)
}
allele_dist <- get_allele_distribtion_dataframe(with_other = TRUE)
allele_dist_all <- get_allele_distribtion_dataframe(FALSE)

# Return the alleles in order of the frequency, except with "Other" always
# at the end.
get_factor_levels <- function(allele, freq) {
    lvls <- allele[order(-freq)]
    if (any(allele == "Other")) {
        lvls <- c(lvls[lvls != "Other"], "Other")
    }
    return(lvls)
}

# Make a bar plot for the distribution of the alleles.
# Provide a value to `max_freq` to set the y-axis limit.
make_allele_dist_barplot <- function(cancer, data,
                                     max_freq = NA) {

    data <- data %>% filter(ras_allele != "WT")

    factor_levels <- get_factor_levels(data$ras_allele, data$allele_freq)

    p <- data %>%
        mutate(ras_allele = factor(ras_allele, levels = !!factor_levels)) %>%
        ggplot(aes(x = ras_allele, y = allele_freq)) +
        geom_col(aes(fill = ras_allele)) +
        scale_fill_manual(
            values = short_allele_pal,
            guide = FALSE) +
        scale_y_continuous(
            limits = c(0, max_freq),
            expand = expand_scale(mult = c(0, 0.02)),
            breaks = round_breaks()
        ) +
        theme_bw(base_size = 8, base_family = "Arial") +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.title = element_blank()
        ) +
        coord_flip() +
        labs(
            title = cancer,
            y = "frequency"
        )

    return(p)
}

# Make a stacked bar plot of the allele frequency.
make_allele_stackedplot <- function(cancer, data, ...) {
    factor_levels <- get_factor_levels(data$ras_allele, data$allele_freq)
    factor_levels <- c("WT", factor_levels[factor_levels != "WT"])
    p <- data %>%
        mutate(ras_allele = factor(ras_allele, levels = !!factor_levels)) %>%
        ggplot(aes(x = "KRAS", y = allele_freq)) +
        geom_col(aes(fill = ras_allele), position = "stack") +
        scale_fill_manual(
            values = short_allele_pal,
            guide = FALSE) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_bw(base_size = 8, base_family = "Arial") +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank()
        )

    return(p)
}

# A wrapper to save the bar plot for a cancer.
save_allele_dist_barplot <- function(barplot, cancer,
                                     size = NA, width = NA, height = NA,
                                     ...) {
    ggsave_wrapper(
        barplot,
        plot_path(GRAPHS_DIR,
                  glue("allele_dist_barplot_{cancer}.svg")),
        size, width, height
    )
}


max_freq <- allele_dist %>%
    filter(ras_allele != "WT") %>%
    pull(allele_freq) %>%
    max()

plots <- allele_dist %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(
        barplot = purrr::map2(cancer, data, make_allele_dist_barplot,
                              max_freq = NA),
        stackedplot = purrr::map2(cancer, data, make_allele_stackedplot)
    ) %>%
    pwalk(save_allele_dist_barplot, width = 4, height = 2.5)

ggsave_wrapper(
    wrap_plots(plots$barplot),
    plot_path(GRAPHS_DIR,
              glue("allele_dist_barplot_all.svg")),
    width = 8, height = 4.5
)



plot_distribution_and_stacked <- function(cancer, data, with_extra_space = FALSE) {
    bp <- make_allele_dist_barplot(cancer, data)
    sp <- make_allele_stackedplot(cancer, data)

    p <- bp + sp + plot_layout(widths = c(20, 1))

    if (with_extra_space) {
        p <- p + plot_spacer() + plot_layout(widths = c(20, 1, 1))
    }

    return(p)
}

plots <- allele_dist %>%
    group_by(cancer) %>%
    nest() %>%
    ungroup() %>%
    arrange(cancer) %>%
    mutate(with_extra_space = rep(c(TRUE, FALSE), 2)) %>%
    pmap(plot_distribution_and_stacked)
ggsave_wrapper(
    patchwork::wrap_plots(plots),
    plot_path(GRAPHS_DIR,
              glue("allele_dist_barplot_stackplot.svg")),
    width = 8, height = 4.5
)

saveRDS(
    plots,
    get_fig_proto_path("allele_dist_barplot_stackplot", 1),
)


# Plot all of the alleles (with at least 3 appearences)
# Need to add the new colors to the palette - I just changed the
#   `short_allele_pal`, though, if this analysis gets any bigger, then I will
#   parameterize the palette for the plotting functions (or pass in `...`).
{
    ORIGINAL_PAL <- short_allele_pal
    added_alleles <- setdiff(unique(allele_dist_all$ras_allele), names(short_allele_pal))
    short_allele_pal <- c(
        short_allele_pal,
        rep(short_allele_pal["Other"], length(added_alleles))
    )
    names(short_allele_pal) <- c(names(ORIGINAL_PAL), added_alleles)

    plots <- allele_dist_all %>%
        filter(num_allele_samples > 2) %>%
        group_by(cancer) %>%
        nest() %>%
        ungroup() %>%
        arrange(cancer) %>%
        mutate(with_extra_space = rep(c(TRUE, FALSE), 2)) %>%
        pmap(plot_distribution_and_stacked)
    ggsave_wrapper(
        patchwork::wrap_plots(plots),
        plot_path(GRAPHS_DIR,
                  glue("allele_dist_barplot_stackplot_all.svg")),
        width = 8, height = 4.5
    )

    saveRDS(
        plots,
        get_fig_proto_path("allele_dist_barplot_stackplot_all", 1, supp = TRUE)
    )

    short_allele_pal <- ORIGINAL_PAL
}


#### ---- Lollipop plot of KRAS mutations ---- ####

codons_to_label <- c(12, 13, 61, 146)

# 'maftools' lollipop plot.
kras_maf <- cancer_full_coding_muts_maf %>%
    filter(hugo_symbol == "KRAS") %>%
    maftools::read.maf(verbose = FALSE)

svg(plot_path(GRAPHS_DIR, "lollipop-kras.svg"), width = 6, height = 4)
maftools::lollipopPlot(kras_maf,
                       "KRAS",
                       AACol = "amino_position",
                       labelPos = codons_to_label,
                       titleSize = c(0.1, 0.1),
                       pointSize = 1.2,
                       axisTextSize = c(0.75, 0.75),
                       showDomainLabel = FALSE)
dev.off()



# My lollipop plot.
kras_lollipop_plot <- cancer_full_coding_muts_maf %>%
    filter(!is_hypermutant & hugo_symbol == "KRAS") %>%
    mutate(amino_position = as.numeric(amino_position)) %>%
    filter(!is.na(amino_position)) %>%
    group_by(cancer, amino_position) %>%
    summarise(num_amino_position = n_distinct(tumor_sample_barcode)) %>%
    ungroup() %>%
    mutate(log_amino_position = log(num_amino_position)) %>%
    group_by(amino_position) %>%
    mutate(
        total_num_amino_position = sum(num_amino_position),
        total_log_amino_position = sum(log_amino_position),
        point_label = ifelse(
            amino_position %in% !!codons_to_label,
            as.character(amino_position),
            NA
        )
    ) %>%
    ungroup() %>%
    ggplot(aes(x = amino_position)) +
    geom_col(
        aes(y = num_amino_position,
            fill = cancer)
    ) +
    geom_point(
        aes(y = total_num_amino_position,
            color = log10(total_num_amino_position)),
        size = 1
    ) +
    geom_text(
        aes(label = point_label,
            y = total_num_amino_position),
        family = "Arial",
        size = 2,
        hjust = 0,
        nudge_x = 5,
        nudge_y = 2
    ) +
    scale_fill_manual(
        values = cancer_palette,
        guide = guide_legend(
            title = NULL,
            label.hjust = 0,
            keywidth = unit(2, "mm"),
            keyheight = unit(2, "mm"),
            ncol = 1
        )
    ) +
    scale_color_viridis_c(
        begin = 0.3, end = 0.9,
        option = "A",
        guide = FALSE
    ) +
    scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, 6000),
        breaks = c(0, 10, 100, 500, 1000, 2000, 4000, 6000)
    ) +
    theme_bw(base_size = 8, base_family = "Arial") +
    theme(
        legend.position = c(0.8, 0.8),
        legend.spacing.x = unit(1, "mm"),
        plot.margin = unit(c(1, 1, 1, 1), "mm")
    ) +
    coord_trans(y = my_trans_log10) +
    labs(
        x = "KRas amino acid sequence",
        y = "log10( count )"
    )

ggsave_wrapper(
    kras_lollipop_plot,
    plot_path(GRAPHS_DIR, "lollipop-kras_2.svg"),
    width = 3, height = 2
)

# Save for use in Figure 1.
saveRDS(
    kras_lollipop_plot,
    get_fig_proto_path("lollipop-kras_2", 1)
)


#### ---- Table of the distribution of alleles ---- ####

alleles_df %<>%
    filter(cancer != "SKCM" & !is_hypermutant) %>%
    mutate(kras_allele = str_remove(ras_allele, "KRAS_"),
           codon = str_extract(ras_allele, "[:digit:]+|WT"))

# Calculate the frequency of the alleles of KRAS by cancer.
calc_frequency_of_alleles_by_cancer <- function(df) {
    df %>%
        group_by(cancer) %>%
        mutate(
            num_cancer_samples = n_distinct(tumor_sample_barcode)
        ) %>%
        group_by(cancer, ras, kras_allele, num_cancer_samples) %>%
        summarise(num_allele_samples = n_distinct(tumor_sample_barcode)) %>%
        ungroup() %>%
        mutate(allele_frequency = num_allele_samples / num_cancer_samples) %>%
        arrange(cancer, -allele_frequency) %>%
        select(cancer, kras_allele,
               num_allele_samples, num_cancer_samples, allele_frequency)
}

# Frequency of each allele across cancers.
alleles_df %>%
    calc_frequency_of_alleles_by_cancer() %>%
    write_tsv(table_path(GRAPHS_DIR, "kras-allele-distribution.tsv"))

# Frequency of each allele across cancers without WT.
alleles_df %>%
    filter(kras_allele != "WT") %>%
    calc_frequency_of_alleles_by_cancer() %>%
    write_tsv(table_path(GRAPHS_DIR, "kras-allele-distribution-noWT.tsv"))


# Frequency of each allele for all cancers without WT.
alleles_df %>%
    mutate(cancer = "ALL") %>%
    calc_frequency_of_alleles_by_cancer() %>%
    write_tsv(table_path(GRAPHS_DIR, "kras-allele-distribution-all.tsv"))

# Frequency of each allele for all cancers without WT.
alleles_df %>%
    mutate(cancer = "ALL") %>%
    filter(kras_allele != "WT") %>%
    calc_frequency_of_alleles_by_cancer() %>%
    write_tsv(table_path(GRAPHS_DIR, "kras-allele-distribution-all-noWT.tsv"))





#' Calculate the frequency of mutations of KRAS codons by cancer.
calc_frequency_of_codons_by_cancer <- function(df) {
    df %>%
        group_by(cancer) %>%
        mutate(
            num_cancer_samples = n_distinct(tumor_sample_barcode)
        ) %>%
        group_by(cancer, codon, num_cancer_samples) %>%
        summarise(num_codon_samples = n_distinct(tumor_sample_barcode)) %>%
        ungroup() %>%
        mutate(codon_frequency = num_codon_samples / num_cancer_samples) %>%
        arrange(cancer, -codon_frequency) %>%
        select(cancer, codon,
               num_codon_samples, num_cancer_samples, codon_frequency)
}

# The frequency of each codon across cancers.
alleles_df %>%
    calc_frequency_of_codons_by_cancer() %>%
    write_tsv(table_path(GRAPHS_DIR, "kras-codon-distribution.tsv"))

# The frequency of each codon across cancers without WT.
alleles_df %>%
    filter(kras_allele != "WT") %>%
    calc_frequency_of_codons_by_cancer() %>%
    write_tsv(table_path(GRAPHS_DIR, "kras-codon-distribution-noWT.tsv"))


# The frequency of codon with all cancers combined.
alleles_df %>%
    mutate(cancer = "ALL") %>%
    calc_frequency_of_codons_by_cancer() %>%
    write_tsv(table_path(GRAPHS_DIR, "kras-codon-distribution-all.tsv"))

# The frequency of codon with all cancers combined without WT.
alleles_df %>%
    filter(kras_allele != "WT") %>%
    mutate(cancer = "ALL") %>%
    calc_frequency_of_codons_by_cancer() %>%
    write_tsv(table_path(GRAPHS_DIR, "kras-codon-distribution-all-noWT.tsv"))
