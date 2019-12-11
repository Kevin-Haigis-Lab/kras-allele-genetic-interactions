
#' A 'ggplot2' based replacement for some of the maftools plots.


#'
#' LAYOUT DIAGRAM FOR `ggoncoplot`:
#'
#'                   o
#'        +---------------------+
#'        |                     |
#'      m |         A           |
#'        |                     |
#'        +---------------------+     p
#'        +---------------------+  +-----+
#'        |                     |  |     |
#'        |                     |  |     |
#'        |                     |  |     |
#'      n |          B          |  |  C  |
#'        |                     |  |     |
#'        |                     |  |     |
#'        |                     |  |     |
#'        |                     |  |     |
#'        +---------------------+  +-----+
#'
#'        +---------------------+
#'      x |          D          |
#'        +---------------------+
#'
#'        +------------------------------+
#'      q |          guide_area          |
#'        +------------------------------+
#'



#' A 'ggplot2' based implementation of 'maftools::oncoplot'.
ggoncoplot <- function(maf,
                       genes = NULL,
                       sampleOrder = NULL,
                       top = 20,
                       cohortSize = NULL,
                       keepGeneOrder = TRUE,
                       removeNonMutated = TRUE,
                       clinicalFeatures = NULL,
                       annotation_pal = NULL,
                       patchwork_layout_params = NULL) {
    # Extract the various data frames from a MAF
    mutations_df <- as_tibble(maf@data)
    variants_per_sample <- as_tibble(maf@variants.per.sample)
    variants_per_gene <- as_tibble(maf@gene.summary)
    clinical_df <- as_tibble(maf@clinical.data)

    # Get cohort size
    if (!is.null(cohortSize)) {
        cohort_size <- cohortSize
    } else {
        cohort_size <- n_distinct(variants_per_sample$Tumor_Sample_Barcode)
    }

    gene_order <- unique(genes)
    if (!keepGeneOrder) {
        gene_order <- sort_genes(maf)
    }

    if (!is.null(top) & !is.infinite(top)) {
        gene_order <- gene_order[seq(1, top)]
        genes <- genes[genes %in% gene_order]
    }

    sample_order <- sampleOrder
    if (is.null(sample_order)) {
        sample_order <- sort_samples(maf, gene_order)
    }

    mutations_df %<>% filter(Hugo_Symbol %in% !!genes)
    variants_per_gene %<>% filter(Hugo_Symbol %in% !!genes)

    variants_per_gene_long <- variants_per_gene %>%
        select(-total, -MutatedSamples, -AlteredSamples) %>%
        pivot_longer(-Hugo_Symbol,
                     names_to = "Variant_Classification",
                     values_to = "num")

    x <- 1
    if (is.null(patchwork_layout_params)) {
        m <- 5 * x
        n <- 20 * x
        o <- 25 * x
        p <- 5 * x
        q <- 8 * x
    } else {
        pm <- patchwork_layout_params
        m <- pm[["m"]] * x
        n <- pm[["n"]] * x
        o <- pm[["o"]] * x
        p <- pm[["p"]] * x
        q <- pm[["q"]] * x
        rm(pm)
    }


    patch_layout <- c(
        area(t = 1, l = 1, b = m, r = o),
        area(t = m+1, l = 1, b = m+n+1, r = o),
        area(t = m+1, l = o+1, b = m+n+1, r = p+o+1),
        area(t = m+n+2, l = 1, b = m+n+x+2, r = o),
        area(t = m+n+x+3, l = 1, b = m+n+x+q+3, r = p+o+1)
    )

    A <- make_top_bar_plot(variants_per_sample,
                           sample_order)
    B <- make_main_tile_plot(mutations_df,
                             gene_order,
                             sample_order)
    C <- make_right_bar_plot(variants_per_gene_long,
                             maf,
                             cohort_size,
                             gene_order)
    if (!is.null(clinicalFeatures)) {
        D <- make_clinical_feature_strip(clinical_df,
                                         sample_order,
                                         annotation_pal)
    } else {
        D <- plot_spacer()
    }


    plot_widths <- c(3, 1)

    onco_patch <- A + B + C + D + guide_area() +
        plot_layout(design = patch_layout,
                    guides = "collect")
    # onco_patch <- (A + plot_spacer() + plot_layout(widths = plot_widths)) /
    #     (B + C + plot_layout(widths = plot_widths)) /
    #     (D + plot_spacer() + plot_layout(widths = plot_widths)) /
    #     guide_area() +
    #     plot_layout(guides = "collect", heights = c(1, 2, 0.1, 1))

    return(onco_patch)
}

####### DELETE ME!
# {
#     ggonco_dir <- "lib/ggmaftools_test_objs"
#     p <- ggoncoplot(maf = readRDS(file.path(ggonco_dir, "maf.rds")),
#                     genes = readRDS(file.path(ggonco_dir, "genes.rds")),
#                     sampleOrder = NULL,
#                     top = 20,
#                     cohortSize = NULL,
#                     keepGeneOrder = TRUE,
#                     removeNonMutated = TRUE,
#                     clinicalFeatures = NULL,
#                     annotation_pal = readRDS(file.path(ggonco_dir, "annotation_pal.rds")))
#     ggsave(
#         file.path(save_dir, "run_ggoncoplot.svg"),
#         p,
#         width = 10, height = 7, unit = "cm"
#     )
# }
#######


#### ---- Plotting components ---- ####

#' Main heat/tile plot of mutations for each gene (y) in each sample (x).
make_main_tile_plot <- function(mutations_df, gene_order, sample_order) {
    p <- mutations_df %>%
        mutate(
            Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode,
                                          levels = sample_order),
            Hugo_Symbol = factor(Hugo_Symbol, levels = rev(gene_order)),
            Variant_Classification = as.character(Variant_Classification),
            Variant_Classification = format_var_names(Variant_Classification)
        ) %>%
        ggplot(aes(x = Tumor_Sample_Barcode, y = Hugo_Symbol)) +
        geom_tile(
            aes(fill = Variant_Classification),
            color = NA
        ) +
        scale_fill_manual(
            values = mod_variant_pal,
            guide = guide_legend(title = NULL, nrow = 1)
        ) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        theme_bw(
            base_size = 6,
            base_family = "Arial"
        ) +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(hjust = 1),
            axis.title = element_blank(),
            legend.position = "bottom",
            axis.ticks = element_blank(),
            legend.key.size = unit(2, "mm"),
            legend.text = element_text(size = 5),
            plot.margin = margin(0, 0, 0, 0, "mm"),
            legend.margin = margin(0, 0, -3, 0, "mm")
        )
    return(p)
}


#' A bar plot along the top with the number of mutations per sample.
make_top_bar_plot <- function(variants_per_sample, sample_order) {
    p <- variants_per_sample %>%
        mutate(
            Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode,
                                          levels = sample_order)
        ) %>%
        ggplot(aes(x = Tumor_Sample_Barcode, y = Variants)) +
        geom_col(fill = "grey30") +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(
            expand = expand_scale(mult = c(0, 0.02)),
            breaks = integer_breaks(rm_vals = c(0))
        ) +
        theme_classic(
            base_size = 6,
            base_family = "Arial"
        ) +
        theme(
            axis.text.x = element_blank(),
            axis.text.y = element_text(hjust = 1),
            axis.title = element_blank(),
            legend.position = "none",
            axis.ticks.x = element_blank(),
            panel.background = element_rect(fill = NA, color = NA),
            panel.grid.major.y = element_line(color = "grey90",
                                              linetype = 2,
                                              size = 0.2),
            plot.margin = margin(0, 0, 0, 0, "mm")
        )
    return(p)

}


#' A bar plot along the right with the number of mutations per gene.
make_right_bar_plot <- function(variants_per_gene_long, maf, cohort_size, gene_order) {
    total_variants_per_gene <- as_tibble(maf@gene.summary) %>%
        select(Hugo_Symbol, total) %>%
        filter(Hugo_Symbol %in% !!gene_order) %>%
        mutate(
            percent_mut = round(total / !!cohort_size * 100),
            percent_mut = paste0(percent_mut, "%"),
            arrange_order = map_int(Hugo_Symbol, ~ which(gene_order == .x))
        ) %>%
        arrange(arrange_order)

    p <- variants_per_gene_long %>%
        mutate(
            Hugo_Symbol = factor(Hugo_Symbol, levels = rev(gene_order)),
            Variant_Classification = format_var_names(Variant_Classification)
        ) %>%
        ggplot(aes(x = Hugo_Symbol, y = num)) +
        geom_col(
            aes(fill = Variant_Classification),
            position = "stack",
            width = 0.95
        ) +
        geom_text(
            data = total_variants_per_gene,
            mapping = aes(x = Hugo_Symbol, y = total, label = percent_mut),
            size = 1.5,
            family = "Arial",
            hjust = 0,
            nudge_y = max(total_variants_per_gene$total) * 0.02
        ) +
        scale_fill_manual(
            values = mod_variant_pal,
            guide = FALSE
        ) +
        scale_x_discrete(
            expand = c(0, 0)
        ) +
        scale_y_continuous(
            expand = expand_scale(mult = c(0, 0.02)),
            limits = c(0, max(total_variants_per_gene$total) * 1.2),
            breaks = integer_breaks(rm_vals = c(0))
        ) +
        theme_classic(
            base_size = 6,
            base_family = "Arial"
        ) +
        theme(
            axis.ticks.y = element_blank(),
            axis.title = element_blank(),
            axis.text.y = element_blank(),
            panel.grid.major.x = element_line(color = "grey90",
                                              linetype = 2,
                                              size = 0.2),
            plot.margin = margin(0, 0, 0, 0, "mm")
        ) +
        coord_flip()
    return(p)
}


#' A strip below to show a clinical feature.
make_clinical_feature_strip <- function(clinical_df, sample_order, pal) {
    colnames(clinical_df)[[2]] <- "a"
    clinical_df$a <- factor(clinical_df$a, levels = names(pal))

    p <- clinical_df %>%
        mutate(
            Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode,
                                          levels = sample_order)
        ) %>%
        ggplot(aes(x = Tumor_Sample_Barcode, y = "")) +
        geom_tile(
            aes(fill = a),
            color = NA
        ) +
        scale_fill_manual(
            values = pal,
            guide = guide_legend(
                title = NULL,
                nrow = 2,
                label.hjust = 0
            )
        ) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        theme_void(
            base_size = 6,
            base_family = "Arial"
        ) +
        theme(
            legend.key.size = unit(2, "mm"),
            legend.text = element_text(size = 5),
            plot.margin = margin(0, 0, 0, 0, "mm"),
            legend.margin = margin(0, 0, -4, 0, "mm")
        ) +
        labs(
            y = "tumor samples"
        )
    return(p)
}



#### ---- Names and colors of variants ---- ####

#' Format the `Variant_Classification` terms to lower case and no "_".
format_var_names <- function(x) {
    y <- str_to_lower(str_replace_all(x, "_", " "))
    y <- str_remove_all(y, " mutation")
    y <- str_replace(y, "in frame", "in-frame")
    y <- str_replace(y, "ins", "ins.")
    y <- str_replace(y, "del", "del.")
    return(y)
}


#' Color palette for variants from 'mafools'.
#' https://github.com/PoisonAlien/maftools/blob/50318fd437fe299244f0f919659eeb8c6246d70f/R/oncomatrix.R
# maftools_variant_pal <- c(
#     RColorBrewer::brewer.pal(12, name = "Paired"),
#     RColorBrewer::brewer.pal(11, name = "Spectral")[1:3],
#     'black', 'violet', 'royalblue', '#7b7060'
# )

maftools_variant_pal <- c(
    "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
    "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928",
    "#9E0142", "#D53E4F", "#F46D43", "black", "violet", "royalblue", "#7b7060"
)

names(maftools_variant_pal) <-c(
    'Nonstop_Mutation', 'Frame_Shift_Del', 'IGR', 'Missense_Mutation','Silent',
    'Nonsense_Mutation', 'RNA','Splice_Site', 'Intron', 'Frame_Shift_Ins',
    'Nonstop_Mutation', 'In_Frame_Del', 'ITD', 'In_Frame_Ins',
    'Translation_Start_Site', "Multi_Hit", 'Amp', 'Del', 'Complex_Event'
)

mod_variant_pal <- maftools_variant_pal
names(mod_variant_pal) <- format_var_names(names(mod_variant_pal))


#### ---- Sorting subroutines ---- ####

#' Sort genes by frequency of mutation.
sort_genes <- function(maf) {
    df <- maf@gene.summary
    return(df$Hugo_Symbol[order(df$total, decreasing = TRUE)])
}


#' If there is a KRAS allele in `genes` return it, otherwise return "WT".
any_kras_allele <- function(genes) {
    kras_gene <- genes[str_detect(genes, "KRAS ")]
    if (length(kras_gene) > 0) { return(kras_gene[[1]]) }
    else { return("WT") }
}


#' Sort samples using a weighted metric.
#' The first weight is by KRAS allele, with "KRAS_other" at the end.
#' The second is by the order of the genes in `ordered_genes`.
sort_samples <- function(maf, ordered_genes) {
    gene_score_tib <- tibble(Hugo_Symbol = rev(ordered_genes)) %>%
        mutate(value = 2 ^ ((1:n()) - 1))

    kras_tib <- tibble(
        kras_allele = rev(ordered_genes[str_detect(ordered_genes, "KRAS")])
    ) %>%
        mutate(kras_value = 1:n())

    maf@data %>%
        as_tibble() %>%
        select(Hugo_Symbol, Tumor_Sample_Barcode) %>%
        left_join(gene_score_tib, by = "Hugo_Symbol") %>%
        group_by(Tumor_Sample_Barcode) %>%
        summarise(score = sum(value),
                  kras_allele = any_kras_allele(Hugo_Symbol)) %>%
        ungroup() %>%
        left_join(kras_tib, by = "kras_allele") %>%
        mutate(kras_value = ifelse(is.na(kras_value), 0, kras_value)) %>%
        arrange(-kras_value, -score) %>%
        pull(Tumor_Sample_Barcode)
}
