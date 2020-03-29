# Compare the subnetworks from COAD alleles.

GRAPHS_DIR <- "40_15_comparing-COAD-allele-subnetworks"
reset_graph_directory(GRAPHS_DIR)

make_depmap_gene_clusters_pairwise_df()
prepare_simple_combined_ppi_gr()


#### ---- Comparing of COAD allele overlap networks ---- ####

annotate_merged_graph <- function(gr, df) {
    for (i in 1:nrow(df)) {
        allele_i <- df$allele[[i]]
        nodes <- igraph::V(df$ppi[[i]])$name
        gr <- gr %N>%
            mutate(allele = case_when(
                name %in% !!nodes & allele == "" ~ !!allele_i,
                name %in% !!nodes ~ paste(allele, !!allele_i, sep = ", "),
                TRUE ~ allele
            ))
    }
    return(gr)
}


make_gray_pal <- function(max_ct, start = 0.8, end = 0.5) {
    pal <- grDevices::gray.colors(max_ct - 1, start = start, end = end)
    names(pal) <- paste(seq(2, max_ct), "alleles", sep = " ")
    return(pal)
}


make_brown_pal <- function(max_ct) {
    pal <- colorRampPalette(c("burlywood2", "burlywood4"))(max_ct)
    names(pal) <- paste(seq(2, max_ct), "alleles", sep = " ")
    return(pal)
}


factor_multiple_node_colors <- function(nc) {
    alleles <- unique(nc[!str_detect(nc, "allele")])
    grays <- unique(nc[str_detect(nc, "allele")])
    levels <- c(sort(alleles), sort(grays))
    return(factor(nc, levels = levels))
}


get_allele_count_label <- function(allele) {
    paste(str_count(allele, ",") + 1, "alleles", sep = " ")
}


annotate_edges_with_clustering <- function(gr) {
    mod_gr <- gr %N>%
        mutate(group = group_label_prop()) %E>%
        mutate(edge_color = ifelse(
            .N()$node_color[to] == .N()$node_color[from],
            as.character(.N()$node_color[to]), NA
        ))
    return(mod_gr)
}


add_group_annotations <- function(gr, anno_tib) {
    mod_gr <- gr %N>%
        mutate(grp_name = NA,
               grp_fill = NA)

    for (i in seq(1, nrow(anno_tib))) {
        nodes <- unlist(anno_tib[i, ]$nodes)
        gn <- anno_tib[i, ]$name[[1]]
        gf <- anno_tib[i, ]$fill[[1]]
        mod_gr <- mod_gr %N>%
            mutate(grp_name = ifelse(name %in% !!nodes, !!gn, grp_name),
                   grp_fill = ifelse(name %in% !!nodes, !!gf, grp_fill))
    }
    return(mod_gr)
}


add_ggforce_annotations <- function(p, anno_tib) {
    anno_pal <- anno_tib %>%
        select(name, fill) %>%
        dplyr::rename(value = fill) %>%
        deframe()
    p <- p +
        ggforce::geom_mark_hull(
            aes(x, y,
                fill = grp_name,
                label = grp_name,
                filter = !is.na(grp_name)),
            color = NA,
            alpha = 0.2,
            label.family = "arial",
            label.fontsize = 7,
            con.cap = unit(1, "mm"),
            con.type = "straight",
            label.buffer = unit(3, "mm"),
            label.fill = NULL,
            concavity = 20,
            expand = unit(3, "mm")
        ) +
        scale_fill_manual(
            values = anno_pal,
            guide = FALSE
        )
    return(p)
}


prepare_overlap_comparison_graph_plot <- function(gr,
                                                  special_labels = NULL,
                                                  annotation_tib = NULL) {
    mod_gr <- gr %>%
        mutate(
            node_color = case_when(
                name %in% !!special_labels ~ "special",
                str_detect(allele, ",") ~ get_allele_count_label(allele),
                TRUE ~ allele
            ),
            node_color = factor_multiple_node_colors(node_color)
        ) %>%
        annotate_edges_with_clustering() %E>%
        mutate(edge_color = factor_multiple_node_colors(edge_color))

    if (!is.null(annotation_tib)) {
        mod_gr <- add_group_annotations(mod_gr, annotation_tib)
    }

    return(mod_gr)
}


make_node_and_edge_palettes <- function(gr) {
    max_ct <- max(str_count(igraph::V(gr)$allele, ",")) + 1
    pal <- c(
        short_allele_pal[names(short_allele_pal) %in% igraph::V(gr)$allele],
        "special" = "indianred1"
    )

    if (max_ct >= 2) {
        pal <- c(pal, make_brown_pal(max_ct))
    }

    edge_pal <- pal[names(pal) %in% levels(igraph::E(gr)$edge_color)]

    return(list(
        node_pal = pal,
        edge_pal = edge_pal
    ))
}




plot_overlap_comparison_graph <- function(gr,
                                           annotation_tib = NULL,
                                           graph_layout = "kk",
                                           node_label_size = 2,
                                           node_label_repel = TRUE,
                                           node_size = 3
                                       ) {

    pals <- make_node_and_edge_palettes(gr)

    p <- ggraph(gr, layout = graph_layout)

    if (!is.null(annotation_tib)) {
        p <- add_ggforce_annotations(p, annotation_tib)
    }

    p <- p +
        geom_edge_link(
            aes(color = edge_color),
            width = 1, alpha = 0.5
        ) +
        geom_node_point(
            aes(color = node_color),
            size = node_size
        ) +
        geom_node_text(
            aes(label = name),
            size = node_label_size,
            family = "Arial",
            repel = node_label_repel
        ) +
        scale_color_manual(
            values = pals$node_pal,
            guide = guide_legend(
                ncol = 2,
                title.hjust = 0.5,
                label.hjust = 0
            )
        ) +
        scale_edge_color_manual(
            values = pals$edge_pal,
            na.value = "grey70",
            guide = FALSE
        ) +
        theme_graph() +
        theme(legend.position = "bottom")
    return(p)
}


print_node_names <- function(gr) {
    # cat(igraph::V(gr)$name, sep = "\n")
    return(gr)
}


extract_common_graph_from_ppi <- function(grs) {
    nodes <- unlist(purrr::map(grs, ~ igraph::V(.x)$name))
    gr <- simple_combined_ppi_gr %N>%
        filter(name %in% !!nodes)
    return(gr)
}


make_overlap_comparison_graph <- function(df) {
    merged_gr <- extract_common_graph_from_ppi(df$ppi) %>%
        convert(to_simple, .clean = TRUE) %>%
        convert(to_undirected, .clean = TRUE) %N>%
        select(name) %E>%
        select(from, to) %N>%
        mutate(allele = "") %>%
        annotate_merged_graph(df)
    return(merged_gr)
}


COAD_GRAPH_ANNOTATIONS <- list(
        "adhesion" = c("MYH6", "MYH9", "L1CAM", "LAMA1", "LAMA2",
                       "DMD", "ITGA6"),
        "extracellular signal\ntransduction" = c("DSCAM", "DCC",
                                                 "ELMO1", "VAV1", "ARHGAP5"),
        "cell morphology" = c("UNK", "AMOT", "USP9X"),
        "immune signaling" = c("NLRP3", "PYCARD", "TOLLIP", "DYNC1H1")
    ) %>%
    enframe(name = "name", value = "nodes") %>%
    mutate(fill = "grey50")


special_nodes <- c("KRAS", "BRAF", "NRAS", "PIK3CA", "APC", "TP53")
genes_to_ignore <- c("TTN")

set.seed(0)
coad_overlap_comparison_plot <- tibble(
        cancer = c("COAD", "COAD", "COAD"),
        allele = c("G12D", "G12V", "G13D")
    ) %>%
    mutate(
        ppi = purrr::map2(cancer, allele,
                          get_overlapped_gr,
                          min_comp_size = 4, ignore_genes = genes_to_ignore),
        data = purrr::map(ppi, ~ .x$data),
        ppi = purrr::map(ppi, ~.x$graph)
    ) %>%
    make_overlap_comparison_graph() %>%
    print_node_names() %>%
    prepare_overlap_comparison_graph_plot(
        special_labels = special_nodes,
        annotation_tib = COAD_GRAPH_ANNOTATIONS
    ) %>%
    plot_overlap_comparison_graph(annotation_tib = COAD_GRAPH_ANNOTATIONS)

ggsave_wrapper(
    coad_overlap_comparison_plot,
    plot_path(GRAPHS_DIR, "coad_overlap_comparison_plot.svg"),
    "large"
)

saveFigRds(coad_overlap_comparison_plot, "coad_overlap_comparison_plot.rds")
