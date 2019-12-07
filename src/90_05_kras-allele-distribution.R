
library(gtable)
library(gridExtra)

library(patchwork)


GRAPHS_DIR <- "90_05_kras-allele-distribution"
reset_graph_directory(GRAPHS_DIR)

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
allele_dist <- alleles_df %>%
    filter(!is_hypermutant) %>%
    mutate(
        ras_allele = str_remove(ras_allele, "KRAS_"),
        ras_allele = fct_other(ras_allele, keep = alleles_to_keep)
    ) %>%
    group_by(cancer) %>%
    mutate(
        ras_allele = as.character(fct_lump(ras_allele, prop = 0.01)),
        num_cancer_samples = n_distinct(tumor_sample_barcode)
    ) %>%
    group_by(cancer, ras, ras_allele, num_cancer_samples) %>%
    summarise(num_allele_samples = n_distinct(tumor_sample_barcode)) %>%
    ungroup()


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

# A temporary data frame to use for plotting. It has a column `allele_freq`
# that has the frequency of each KRAS allele in each cancer.
df <- allele_dist %>%
    filter(cancer != "SKCM") %>%
    mutate(allele_freq = num_allele_samples / num_cancer_samples)

max_freq <- df %>% filter(ras_allele != "WT") %>% pull(allele_freq) %>% max()

plots <- df %>%
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

plots <- df %>%
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




my_trans_log10 <- scales::trans_new(
    name = "log10 pseudo-count +1",
    transform = function(x) { log10(x + 1) },
    inverse = function(x) { exp(x) - 1 }
)


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
        aes(y = num_amino_position, fill = cancer)
    ) +
    geom_point(
        aes(y = total_num_amino_position, color = log10(total_num_amino_position)),
        size = 1
    ) +
    geom_text(
        aes(label = point_label, y = total_num_amino_position),
        family = "Arial",
        hjust = 0,
        nudge_x = 2,
        nudge_y = 1
    ) +
    scale_fill_manual(values = cancer_palette) +
    scale_color_viridis_c(
        begin = 0.3, end = 0.9,
        option = "A"
    ) +
    scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, 6000),
        breaks = c(0, 10, 100, 500, 1000, 2000, 3000, 4000, 6000)
    ) +
    theme_bw(base_size = 8, base_family = "Arial") +
    theme(
        legend.position = "right",
        plot.margin = unit(c(1, 1, 1, 1), "mm")
    ) +
    coord_trans(y = my_trans_log10) +
    labs(
        x = "KRAS amino acid sequence",
        y = "log10( number of mutations + 1 )",
        color = "log10( count )",
        fill = ""
    )
ggsave_wrapper(kras_lollipop_plot,
               plot_path(GRAPHS_DIR, "lollipop-kras_2.svg"),
               width = 5, height = 4)



#### ---- Table of the distribution of alleles ---- ####

# Frequency of each allele across cancers.
alleles_df %>%
    filter(cancer != "SKCM") %>%
    filter(!is_hypermutant) %>%
    mutate(
        kras_allele = str_remove(ras_allele, "KRAS_")
    ) %>%
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
           num_allele_samples, num_cancer_samples, allele_frequency) %>%
    write_tsv(file.path("tables", GRAPHS_DIR, "kras-allele-distribution.tsv"))

# The frequency of each codon across cancers.
alleles_df %>%
    filter(cancer != "SKCM") %>%
    filter(!is_hypermutant) %>%
    mutate(
        codon = str_extract(ras_allele, "[:digit:]+|WT")
    ) %>%
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
           num_codon_samples, num_cancer_samples, codon_frequency) %>%
    write_tsv(file.path("tables", GRAPHS_DIR, "kras-codon-distribution.tsv"))

# The frequency of codon with all cancers combined.
alleles_df %>%
    filter(cancer != "SKCM") %>%
    filter(!is_hypermutant) %>%
    mutate(
        codon = str_extract(ras_allele, "[:digit:]+|WT"),
        num_samples = n_distinct(tumor_sample_barcode)
    ) %>%
    group_by(codon, num_samples) %>%
    summarise(num_codon_samples = n_distinct(tumor_sample_barcode)) %>%
    ungroup() %>%
    mutate(codon_frequency = num_codon_samples / num_samples) %>%
    arrange(-codon_frequency) %>%
    select(codon, num_codon_samples, num_samples, codon_frequency) %>%
    write_tsv(file.path("tables", GRAPHS_DIR,
                        "kras-codon-distribution-allcancers.tsv"))
