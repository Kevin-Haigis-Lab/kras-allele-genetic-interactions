
library(gtable)
library(gridExtra)

alleles_df <- cancer_muts_df %>%
    group_by(tumor_sample_barcode) %>%
    slice(1) %>%
    ungroup() %>%
    select(cancer, tumor_sample_barcode, dataset, target,
           is_hypermutant, ras, ras_allele) %>%
    unique()

alleles_to_keep <- names(short_allele_pal)
alleles_to_keep <- alleles_to_keep[alleles_to_keep != "Other"]

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


get_factor_levels <- function(allele, freq) {
    lvls <- allele[order(-freq)]
    if (any(allele == "Other")) {
        lvls <- c(lvls[lvls != "Other"], "Other")
    }
    return(lvls)
}

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
            expand = expand_scale(mult = c(0, 0.02))
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

save_allele_dist_barplot <- function(barplot, cancer,
                                     size = NA, width = NA, height = NA,
                                     ...) {
    ggsave_wrapper(
        barplot,
        plot_path("90_05_kras-allele-distribution",
                  glue("allele_dist_barplot_{cancer}.svg")),
        size, width, height
    )
}

df <- allele_dist %>%
    filter(cancer != "SKCM") %>%
    mutate(allele_freq = num_allele_samples / num_cancer_samples)

max_freq <- df %>% filter(ras_allele != "WT") %>% pull(allele_freq) %>% max()

plots <- df %>%
    group_by(cancer) %>%
    nest() %>%
    mutate(
        barplot = purrr::map2(cancer, data, make_allele_dist_barplot,
                              max_freq = !!max_freq),
        stackedplot = purrr::map2(cancer, data, make_allele_stackedplot)
    ) %>%
    pwalk(save_allele_dist_barplot, width = 4, height = 2.5)


p <- cowplot::plot_grid(plotlist = plots$barplot, align = "hv", nrow = 2)
cowplot::save_plot(
    plot_path("90_05_kras-allele-distribution",
              glue("allele_dist_barplot_all.svg")),
    plot = p,
    base_width = 8, base_height = 4.5
)

# A place-holder for plotting.
empty_grob <- ggplotGrob(
    model_data  %>%
    sample_n(10) %>%
    ggplot(aes(x = dep_map_id, y = hugo_symbol)) +
    geom_tile(fill = NA, color = NA) +
    theme_void()
)

subplotlist <- rep(list(NA), 4)
names(subplotlist) <- sort(unique(df$cancer))
for (CANCER in unique(df$cancer)) {
    bp <- ggplotGrob(plots$barplot[[which(plots$cancer == CANCER)]])
    sp <- ggplotGrob(plots$stackedplot[[which(plots$cancer == CANCER)]])

    shared_height <- grid::unit.pmax(bp$heights, sp$heights)
    bp$heights <- shared_height
    sp$heights <- shared_height

    subplotlist[[CANCER]] <- arrangeGrob(
        grobs = list(bp, sp),
        nrow = 1,
        widths = c(6, 1)
    )
}

fullplot <- arrangeGrob(
    grobs = subplotlist,
    layout_matrix = rbind(
        c(1, NA, 2),
        c(3, NA, 4)
    ),
    nrow = 2,
    widths = c(10, 1, 10)
    )
ggsave_wrapper(fullplot,
               plot_path("90_05_kras-allele-distribution",
                         glue("allele_dist_barplot_stackplot.svg")),
               "wide")
