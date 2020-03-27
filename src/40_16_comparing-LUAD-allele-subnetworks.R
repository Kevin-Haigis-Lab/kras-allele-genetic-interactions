# Compare the subnetworks from LUAD alleles.

source(file.path("src", "40_15_comparing-COAD-allele-subnetworks.R"))

GRAPHS_DIR <- "40_16_comparing-LUAD-allele-subnetworks"
reset_graph_directory(GRAPHS_DIR)


isolate_kras_subnetwork <- function(gr) {
    gr %N>%
        morph(to_components) %>%
        mutate(has_kras = any(name == "KRAS")) %>%
        unmorph() %>%
        filter(has_kras)
}



aprior_genes <- wide_genetic_interaction_df %>%
    filter(
        hugo_symbol == "KRAS" | !is.na(KEGG) | !is.na(CGC) | !is.na(BioID)
    ) %>%
    pull(hugo_symbol) %>%
    unique()


prep_luad_ggforce <- function(gr) {
    grp1 <- c("NFE2L2", "KEAP1", "DCUN1D3", "ANAPC2")
    gr %N>%
        mutate(ggforce_grp = ifelse(name %in% grp1, "test group", "NA"))
}


add_luad_ggforce <- function(p) {
    p +
        ggforce::geom_mark_hull(
            aes(x, y,
                filter = ggforce_grp == "test group"),
            color = NA,
            fill = "green",
            alpha = 0.15
        )
}



genes_to_ignore <- c("TTN")
special_nodes <- c("KRAS", "BRAF", "NRAS", "PIK3CA", "APC", "TP53",
                   "EGFR", "SMAD4", "SMARCA4", "STK11")

luad_overlap_comparison_plot <- tibble(
        cancer = c("LUAD", "LUAD"),
        allele = c("G12C", "G12V")
    ) %>%
    mutate(
        ppi = map2(cancer, allele,
                   get_overlapped_gr,
                   min_comp_size = 4,
                   ignore_genes = genes_to_ignore),
        data = map(ppi, ~ .x$data),
        ppi = map(ppi, ~.x$graph),
        ppi = map(ppi, isolate_kras_subnetwork),
        ppi = map(ppi, ~ .x %N>% filter(name %in% !!aprior_genes)),
        ppi = map(ppi, get_giant_component)
    ) %>%
    make_overlap_comparison_graph() %>%
    prep_luad_ggforce() %>%
    plot_overlap_comparison_graph2(special_labels = special_nodes) %>%
    add_luad_ggforce()

ggsave_wrapper(
    luad_overlap_comparison_plot,
    plot_path(GRAPHS_DIR, "luad_overlap_comparison_plot.svg"),
    "large"
)

saveFigRds(luad_overlap_comparison_plot, "luad_overlap_comparison_plot")
