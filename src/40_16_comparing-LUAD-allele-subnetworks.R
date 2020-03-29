# Compare the subnetworks from LUAD alleles.

source(file.path("src", "40_15_comparing-COAD-allele-subnetworks.R"))

GRAPHS_DIR <- "40_16_comparing-LUAD-allele-subnetworks"
reset_graph_directory(GRAPHS_DIR)
reset_table_directory(GRAPHS_DIR)


only_allele_clusters <- function(gr, special_nodes = NULL) {
    gr %N>%
        filter(!str_detect(allele, ",") | name %in% !!special_nodes) %E>%
        filter(
            .N()$allele[to] == .N()$allele[from] |
            .N()$allele[to] %in% special_nodes |
            .N()$allele[from] %in% special_nodes
        ) %>%
        activate("nodes")
}



aprior_genes <- wide_genetic_interaction_df %>%
    filter(
        hugo_symbol == "KRAS" | !is.na(KEGG) | !is.na(CGC) | !is.na(BioID)
    ) %>%
    pull(hugo_symbol) %>%
    unique()

LUAD_GRAPH_ANNOTATIONS <- list(
        # G12C
        "WNT signaling" = c("WNT11", "WNT9A", "LRP6", "WNT2", "CTNNB1",
                            "NOTCH1", "FRAT2", "PTPRJ", "TCF7L1", "PYGO1",
                            "CDH10", "JAG2"),
        "ECM" = c("COL5A2", "COL6A6", "COL28A1", "COL24A1", "COL14A1",
                  "COL11A1", "COL3A1", "COL8A2", "FN1", "FLNA", "LAMC1",
                  "VCAN", "ITGA4"),
        "cell cycle\nregulation" = c("MDC1", "MRE11A", "FANCM", "RAD54B",
                                     "EXO1", "SMC4", "NCAPD2", "STAG2", "CENPE",
                                     "SUN1", "AHCTF1"),
        "regulation of ERK1/2" = c("EPHB6", "EPHA5", "EPHA6", "EPHA8", "EPHA2"),
        "chromatin\nremodeling" = c("ARID1B", "ARID1A", "HCFC1", "SRCAP",
                                    "ASXL1"),
        "L1CAM interactions" = c("SCN5A", "SCN10A", "NRCAM", "NFASC", "CNTN1"),
        "glycosylation" = c("ADAMTS16", "ADAMTS8", "ADAMTS15", "THSD7B", "ACAN",
                            "HAPLN1"),
        # G12V
        "vasoconstriction" = c("EDN1", "KALRN", "ADRA1D", "EDNRA", "ADCY2",
                               "EDN1"),
        "transcriptional\nregulation" = c("POLR2A", "HNRNPUL1", "WWP2", "CREB1",
                                          "HDAC1", "ZNF639", "KMT2C", "RUNX1",
                                          "MGA", "ZNF219", "SATB1", "SFMBT1"),
        "metabolism" = c("IDH1", "AK1", "B3GLCT"),
        "ubiquitination" = c("HERC2", "CCNF", "LTN1", "FBXO7"),
        "DSB repair" = c("RIF1", "RAD50", "TERF2IP", "RMI2")

    ) %>%
    enframe(name = "name", value = "nodes") %>%
    mutate(fill = "grey50",
           allele = c(rep("G12C", 7), rep("G12V", n() - 7)))


genes_to_ignore <- c("TTN")
special_nodes <- c("KRAS", "BRAF", "NRAS", "PIK3CA", "APC", "TP53",
                   "EGFR", "SMAD4", "SMARCA4", "STK11")


luad_overlap_comparison_ppin <- tibble(
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
        # ppi = map(ppi, ~ .x %N>% filter(name %in% !!aprior_genes)),
        # ppi = map(ppi, get_giant_component)
    ) %>%
    make_overlap_comparison_graph() %>%
    only_allele_clusters(special_nodes = special_nodes) %>%
    filter_component_size(min_size = 6) %T>%
    print_functional_groups(
        file_name = table_path(GRAPHS_DIR, "enrichr-output.tsv"),
        ignore_genes = special_nodes
    ) %>%
    prepare_overlap_comparison_graph_plot(
        special_labels = special_nodes,
        annotation_tib = LUAD_GRAPH_ANNOTATIONS
    )


for (al in c("G12C", "G12V")) {

    anno_tib <- LUAD_GRAPH_ANNOTATIONS %>%
        filter(allele == !!al)

    fxn_grp_name <- table_path(GRAPHS_DIR,
                               as.character(glue("{al}_enrichr-output.tsv")))
    plt_name <- plot_path(
        GRAPHS_DIR, as.character(glue("luad-{al}_overlap_comparison_plot.svg"))
    )

    node_size <- ifelse(al == "G12C", 1, 3)


    allele_ppin <- luad_overlap_comparison_ppin %N>%
        filter(allele == !!al) %>%
        print_functional_groups(
            file_name = fxn_grp_name,
            ignore_genes = special_nodes
        )

    if (al == "G12C") {
        allele_ppin <- allele_ppin %N>%
            mutate(name = ifelse(name %in% unlist(anno_tib$nodes), name, NA))
    }

    allele_ppin %>%
        plot_overlap_comparison_graph(
            annotation_tib = anno_tib,
            graph_layout = ifelse(al == "G12C", "nicely", "stress"),
            node_size = node_size
        ) %T>%
        ggsave_wrapper(plt_name, "large") %>%
        saveFigRds(as.character(glue("luad-{al}_overlap_comparison_plot")))
}
