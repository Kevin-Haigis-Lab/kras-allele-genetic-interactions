# Rainfall plots for the more interesting and strongest genetic interactions.


GRAPHS_DIR <- "20_50_rainfall-plots"
reset_graph_directory(GRAPHS_DIR)

# Data frames and functions for preparing MAFs from cancer data.
source(file.path("src", "20_49_rainfall-plots-subroutines.R"))


#### ---- Oncoplots for interactors of KRAS by allele ---- ####
# Plot the `maftools::oncostrip()` for the genes mutually exclusive with one
# allele for all of the alleles.


# Get all KRAS alleles for a cancer
get_alleles_in_cancer <- function(cancer) {
  cancer_coding_muts_df %>%
    filter(cancer == !!cancer) %>%
    u_pull(ras_allele)
}


# Sort genes into: KRAS, "KRAS other", the rest by mutation frequency.
sort_genes_elevate_kras <- function(maf) {
  gene_order <- sort_genes(maf)

  kras_allele <- gene_order[str_detect(gene_order, "KRAS ")]
  if (any(kras_allele == "KRAS other")) {
    kras_allele <- c(
      kras_allele[kras_allele != "KRAS other"],
      "KRAS other"
    )
  }
  gene_order <- c(
    kras_allele,
    gene_order[!(gene_order %in% kras_allele)]
  )
  return(gene_order)
}



# Sort the genes for oncostrip and keep KRAS at the top.
sort_maf <- function(maf) {
  gene_order <- sort_genes_elevate_kras(maf)
  sample_order <- sort_samples(maf, ordered_genes = gene_order)
  return(list(
    gene_order = gene_order,
    sample_order = sample_order
  ))
}


# Tibble of number of samples per cancer to be used for the oncoplot.
cancer_count_tib <- cancer_muts_df %>%
  group_by(cancer, tumor_sample_barcode) %>%
  slice(1) %>%
  ungroup() %>%
  unique() %>%
  count(cancer)

# A thin wrapper around `maftools::oncoplot`.
oncoplot_wrapper <- function(maf,
                             gene_and_sample_orders,
                             save_path,
                             top = 10,
                             cancer = NULL,
                             annotate_kras_allele = FALSE,
                             save_size = list(width = 9, height = 5),
                             return_ggonco = FALSE,
                             ggonco_only = FALSE,
                             with_legend_section = TRUE) {
  if (!is.null(cancer)) {
    cohortSize <- unlist(filter(cancer_count_tib, cancer == !!cancer)$n)
  } else {
    cohortSize <- NULL
  }

  if (annotate_kras_allele) {
    idx <- names(short_allele_pal) %in% getClinicalData(maf)$kras_allele
    pal <- short_allele_pal[idx]

    if (!ggonco_only) {
      svg(save_path, width = save_size$width, height = save_size$height)
      oncoplot(
        maf,
        genes = gene_and_sample_orders$gene_order,
        sampleOrder = gene_and_sample_orders$sample_order,
        top = top,
        cohortSize = cohortSize,
        keepGeneOrder = TRUE,
        GeneOrderSort = FALSE,
        removeNonMutated = TRUE,
        clinicalFeatures = "kras_allele",
        annotationColor = list(kras_allele = pal),
        sepwd_samples = 0,
        fontSize = 0.6
      )
      dev.off()
    }

    if (return_ggonco) {
      p <- ggoncoplot(
        maf = maf,
        genes = gene_and_sample_orders$gene_order,
        sampleOrder = gene_and_sample_orders$sample_order,
        top = top,
        cohortSize = cohortSize,
        keepGeneOrder = TRUE,
        removeNonMutated = TRUE,
        clinicalFeatures = "kras_allele",
        annotation_pal = pal,
        with_legend_section = with_legend_section
      )
      return(p)
    }
  } else {
    if (!ggonco_only) {
      svg(save_path, width = 9, height = 5)
      oncoplot(
        maf,
        genes = gene_and_sample_orders$gene_order,
        sampleOrder = gene_and_sample_orders$sample_order,
        top = top,
        cohortSize = cohortSize,
        keepGeneOrder = TRUE,
        GeneOrderSort = FALSE,
        removeNonMutated = TRUE,
        sepwd_samples = 0,
        fontSize = 0.6
      )
      dev.off()
    }

    if (return_ggonco) {
      p <- ggoncoplot(
        maf = maf,
        genes = gene_and_sample_orders$gene_order,
        sampleOrder = gene_and_sample_orders$sample_order,
        top = top,
        cohortSize = cohortSize,
        keepGeneOrder = TRUE,
        removeNonMutated = TRUE,
        annotation_pal = pal,
        with_legend_section = with_legend_section
      )
      return(p)
    }
  }
}


ggoncoplot_wrapper <- function(maf,
                               gene_and_sample_orders,
                               save_path = NULL, # PARAM NOT USED
                               top = 10,
                               cancer = NULL,
                               annotate_kras_allele = FALSE,
                               save_size = NULL) { # PARAM NOT USED

  if (!is.null(cancer)) {
    cohortSize <- unlist(filter(cancer_count_tib, cancer == !!cancer)$n)
  } else {
    cohortSize <- NULL
  }
  if (annotate_kras_allele) {
    idx <- names(short_allele_pal) %in% getClinicalData(maf)$kras_allele
    pal <- short_allele_pal[idx]

    p <- ggoncoplot(
      maf,
      genes = gene_and_sample_orders$gene_order,
      sampleOrder = gene_and_sample_orders$sample_order,
      top = top,
      cohortSize = cohortSize,
      keepGeneOrder = TRUE,
      removeNonMutated = TRUE,
      clinicalFeatures = "kras_allele",
      annotationColor = list(kras_allele = pal)
    )
    return(p)
  } else {
    p <- ggoncoplot(
      maf,
      genes = gene_and_sample_orders$gene_order,
      sampleOrder = gene_and_sample_orders$sample_order,
      top = top,
      cohortSize = cohortSize,
      keepGeneOrder = TRUE,
      removeNonMutated = TRUE
    )
    invisible(p)
  }
}



# MAIN PLOTTING FUNCTION
# Plot the oncoplot for a `cancer` and KRAS `allele` for the
#   top `top_n_genes` genes, ignoring any in `ignore_genes`.
# (The `...` doesn't go anywhere.)
allele_exclusivity_oncoplot <- function(cancer, allele,
                                        save_name_template,
                                        gg_save_name_template = NULL,
                                        interaction_type = "exclusivity",
                                        top_n_genes = 20,
                                        ignore_genes = NULL,
                                        annotate_kras_allele = FALSE,
                                        ...) {
  mutex_genes <- genetic_interaction_gr %N>%
    filter(!name %in% !!ignore_genes) %E>%
    filter(cancer == !!cancer &
      kras_allele == !!allele &
      genetic_interaction %in% !!interaction_type) %>%
    top_n(top_n_genes, -p_val) %N>%
    filter(centrality_degree(mode = "all") > 0) %N>%
    as_tibble() %>%
    filter(!is_kras) %>%
    pull(name)

  if (length(mutex_genes) == 0) {
    return()
  }

  mutex_genes <- c(
    mutex_genes,
    "KRAS"
  )

  maf <- get_maf(mutex_genes, cancer, allele,
    replace_kras_with_allele = TRUE,
    group_other_alleles = TRUE
  )

  gene_and_sample_orders <- sort_maf(maf)

  allele_short <- str_remove(allele, "KRAS_")
  save_path <- plot_path(
    GRAPHS_DIR,
    glue(save_name_template)
  )

  if (!is.null(gg_save_name_template)) {
    gg_save_path <- plot_path(
      GRAPHS_DIR,
      glue(gg_save_name_template)
    )
  } else {
    gg_save_path <- plot_path(
      GRAPHS_DIR,
      paste0("gg_", glue(save_name_template))
    )
  }

  g <- oncoplot_wrapper(
    maf = maf,
    gene_and_sample_orders = gene_and_sample_orders,
    save_path = save_path,
    top = top_n_genes + 2,
    cancer = cancer,
    annotate_kras_allele = annotate_kras_allele,
    return_ggonco = TRUE
  )
  ggsave_wrapper(g, gg_save_path, "wide")
}


# Uninteresting genes
uninteresting_genes <- c("TTN", "FN1")

# Genes to ignore in most situations.
commonly_interacting_genes <- c(
  uninteresting_genes,
  "BRAF", "APC", "TP53", "NRAS", "PIK3CA"
)


cancer_counts_df <- cancer_muts_df %>%
  group_by(cancer, ras_allele, tumor_sample_barcode) %>%
  slice(1) %>%
  ungroup() %>%
  count(cancer, ras_allele) %>%
  dplyr::rename(allele = "ras_allele")

cancer_counts_df %>%
  filter(n > 10) %>%
  pwalk(
    allele_exclusivity_oncoplot,
    save_name_template = "{cancer}_{allele_short}_exclusivity_oncostrip_notCommonInteractors.svg",
    top_n_genes = 15,
    ignore_genes = commonly_interacting_genes,
    annotate_kras_allele = TRUE
  )

cancer_counts_df %>%
  filter(n > 10) %>%
  pwalk(
    allele_exclusivity_oncoplot,
    save_name_template = "{cancer}_{allele_short}_exclusivity_oncostrip_allInteractors.svg",
    top_n_genes = 15,
    ignore_genes = uninteresting_genes,
    annotate_kras_allele = TRUE
  )

cancer_counts_df %>%
  filter(n > 10) %>%
  pwalk(
    allele_exclusivity_oncoplot,
    save_name_template = "{cancer}_{allele_short}_bothInteractions_oncostrip_allInteractors.svg",
    interaction_type = c("comutation", "exclusivity"),
    top_n_genes = 15,
    ignore_genes = uninteresting_genes,
    annotate_kras_allele = TRUE
  )

cancer_counts_df %>%
  filter(n > 10) %>%
  pwalk(
    allele_exclusivity_oncoplot,
    save_name_template = "{cancer}_{allele_short}_comutation_oncostrip_allInteractors.svg",
    interaction_type = "comutation",
    top_n_genes = 15,
    ignore_genes = uninteresting_genes,
    annotate_kras_allele = TRUE
  )



#### ---- Specific genes to  make oncoplots for ---- ####

# MAIN PLOTTING FUNCTION
# Plot the oncoplot for a `cancer` and KRAS `allele` for a specific set of
#   genes. The genes are checked against the interaction graph.
# (The `...` doesn't go anywhere.)
allele_specificgenes_oncoplot <- function(cancer,
                                          kras_allele,
                                          genes,
                                          interaction_type,
                                          save_name,
                                          save_dir,
                                          ignore_genes = NULL,
                                          top_n_genes = 10,
                                          replace_kras_with_allele = TRUE,
                                          annotate_kras_allele = FALSE,
                                          print_missing_genes = TRUE,
                                          keep_gg_onco_proto = FALSE,
                                          ...) {
  genes <- unlist(genes)

  checked_genes <- genetic_interaction_df %>%
    filter(cancer == !!cancer &
      kras_allele == kras_allele &
      hugo_symbol %in% !!genes) %>%
    jhcutils::u_pull(hugo_symbol)

  n_g <- n_distinct(genes)
  n_cg <- n_distinct(checked_genes)
  cat(glue("{cancer}, {kras_allele}: {n_g - n_cg} gene(s) were removed."), "\n")
  if (n_g != n_cg & print_missing_genes) {
    cat(genes[!genes %in% checked_genes], "\n")
  }

  if (length(checked_genes) == 0) {
    cat("**There were zero genes to plot; returning early.\n")
    return()
  }

  checked_genes <- c(
    checked_genes,
    "KRAS"
  )

  maf <- get_maf(checked_genes, cancer, kras_allele,
    replace_kras_with_allele = replace_kras_with_allele,
    group_other_alleles = replace_kras_with_allele
  )

  gene_and_sample_orders <- sort_maf(maf)

  save_path <- plot_path(
    save_dir,
    save_name
  )

  gg_save_path <- plot_path(
    save_dir,
    glue("gg_{save_name}")
  )

  img_height <- ((length(checked_genes) + 1) / 5) + 2

  g <- oncoplot_wrapper(
    maf = maf,
    gene_and_sample_orders = gene_and_sample_orders,
    save_path = save_path,
    top = top_n_genes + 1 + as.numeric(replace_kras_with_allele),
    cancer = cancer,
    annotate_kras_allele = annotate_kras_allele,
    save_size = list(width = 9, height = img_height),
    return_ggonco = TRUE,
    with_legend_section = TRUE
  )

  ggsave_wrapper(g, gg_save_path, "wide")

  g <- oncoplot_wrapper(
    maf = maf,
    gene_and_sample_orders = gene_and_sample_orders,
    save_path = save_path,
    top = top_n_genes + 1 + as.numeric(replace_kras_with_allele),
    cancer = cancer,
    annotate_kras_allele = annotate_kras_allele,
    save_size = list(width = 9, height = img_height),
    return_ggonco = TRUE,
    ggonco_only = TRUE,
    with_legend_section = FALSE
  )

  if (keep_gg_onco_proto) {
    save_gg_onco_proto(g, save_name, cancer)
  }
}


# Save a ggoncoplot (or any plot really) for use in Figure 2.
save_gg_onco_proto <- function(gg_obj, save_name, cancer) {
  n <- file_sans_ext(save_name)
  cat("saving gg object for:", n, "\n")

  if (cancer %in% c("COAD", "LUAD")) {
    saveFigRds(gg_obj, n)
  }

  invisible(NULL)
}


specific_oncoplot_info_tib <- bind_rows(
  tibble(
    cancer = "COAD",
    allele = "A146T",
    interaction_type = "comutation",
    name_suffix = "",
    genes = c("PORCN", "RGS20", "DAPK1"),
  ),
  tibble(
    cancer = "COAD",
    allele = "G12D",
    interaction_type = "exclusivity",
    name_suffix = "",
    genes = c("SCN10A", "DICER1", "SDK1", "TRPM2"),
    keep_gg_onco_proto = TRUE
  ),
  tibble(
    cancer = "COAD",
    allele = "G12D",
    interaction_type = "comutation",
    name_suffix = "",
    genes = c("MAGEC1", "AMER1", "TGIF1"),
    keep_gg_onco_proto = TRUE
  ),
  tibble(
    cancer = "COAD",
    allele = "G12D",
    interaction_type = "exclusiveYWHAZ",
    name_suffix = "",
    genes = c("BRAF", "ZC3H13", "DYNC1H1", "MYH9", "MLLT4"),
    keep_gg_onco_proto = TRUE
  ),
  tibble(
    cancer = "COAD",
    allele = "G12D",
    interaction_type = "comutYWHAZ",
    name_suffix = "",
    genes = c("FAM65B", "IPO7", "CTNNB1", "PCBP1", "PIK3CA", "APC"),
    keep_gg_onco_proto = TRUE
  ),
  tibble(
    cancer = "COAD",
    allele = "G12V",
    interaction_type = "exclusivity",
    name_suffix = "",
    genes = c("DNAH2", "RYR1", "PKHD1", "DIDO1", "PKD1", "KALRN"),
    keep_gg_onco_proto = TRUE
  ),
  tibble(
    cancer = "COAD",
    allele = "G12V",
    interaction_type = "comutation",
    name_suffix = "",
    genes = c("APC", "PIK3CA", "SMAD4", "AMER1", "MCC"),
    keep_gg_onco_proto = TRUE
  ),
  tibble(
    cancer = "COAD",
    allele = "G13D",
    interaction_type = "comutation",
    name_suffix = "",
    genes = c("ANK3", "KIF4B", "SPHKAP")
  ),
  tibble(
    cancer = "COAD",
    allele = "G13D",
    interaction_type = "exclusivity",
    name_suffix = "",
    genes = c("TP53", "CAMTA1", "ZFHX4")
  ),
  tibble(
    cancer = "LUAD",
    allele = "G12C",
    interaction_type = "exclusivity",
    name_suffix = "(MUC)",
    genes = c("MUC4", "MUC17")
  ),
  tibble(
    cancer = "LUAD",
    allele = "G12C",
    interaction_type = "exclusivity",
    name_suffix = "",
    genes = c("CSMD2", "DNAH5", "SMARCA4", "CHD5"),
  ),
  tibble(
    cancer = "LUAD",
    allele = "G12C",
    interaction_type = "comutation",
    name_suffix = "",
    genes = c("PRDM9", "RIMS2", "STK11", "SLITRK2")
  ),
  tibble(
    cancer = "LUAD",
    allele = "G12D",
    interaction_type = "exclusivity",
    name_suffix = "",
    genes = c("TP53")
  ),
  tibble(
    cancer = "LUAD",
    allele = "G12V",
    interaction_type = "exclusivity",
    name_suffix = "",
    genes = c("FAT4", "LAMA2", "SPHKAP", "ZNF804A")
  ),
  tibble(
    cancer = "PAAD",
    allele = "G12D",
    interaction_type = "exclusivity",
    name_suffix = "(RYR)",
    genes = c("RYR3", "RYR2")
  ),
  tibble(
    cancer = "PAAD",
    allele = "G12D",
    interaction_type = "exclusivity",
    name_suffix = "",
    genes = c("RYR3", "RYR2", "TGFBR2", "GNAS")
  ),
  tibble(
    cancer = "PAAD",
    allele = "G12D",
    interaction_type = "comutation",
    name_suffix = "_TP53",
    genes = c("TP53")
  ),
  tibble(
    cancer = "PAAD",
    allele = "G12D",
    interaction_type = "comutation",
    name_suffix = "",
    genes = c("ABCC9", "ZNF831")
  ),
  tibble(
    cancer = "PAAD",
    allele = "G12R",
    interaction_type = "exclusivity",
    name_suffix = "",
    genes = c("RNF43", "DNAH11", "GNAS", "DNAH5", "SMARCA4", "PCDH15")
  ),
  tibble(
    cancer = "PAAD",
    allele = "G12R",
    interaction_type = "comutation",
    name_suffix = "",
    genes = c("BRCA2", "SCN10A")
  ),
  tibble(
    cancer = "PAAD",
    allele = "G12V",
    interaction_type = "comutation",
    name_suffix = "_TP53",
    genes = c("TP53")
  ),
  tibble(
    cancer = "PAAD",
    allele = "G12V",
    interaction_type = "exclusivity",
    name_suffix = "",
    genes = c("GNAS", "FAT4", "HMCN1", "CSMD2")
  ),
  tibble(
    cancer = "PAAD",
    allele = "G12V",
    interaction_type = "comutation",
    name_suffix = "",
    genes = c("RYR3", "RYR2", "PCDHA4")
  ),
  tibble(
    cancer = "PAAD",
    allele = "G12V",
    interaction_type = "exclusivitycomutation",
    name_suffix = "(RYR)",
    genes = c("RYR1", "RYR3", "RYR2")
  ),
  tibble(
    cancer = "PAAD",
    allele = "Q61H",
    interaction_type = "exclusivity",
    name_suffix = "",
    genes = c("SMAD4", "ARID1A", "DNAH5", "ATM", "KDM6A")
  ),
  tibble(
    cancer = "PAAD",
    allele = "Q61H",
    interaction_type = "comutation",
    name_suffix = "",
    genes = c("TP53", "RBM10")
  ),
  tibble(
    cancer = "PAAD",
    allele = "Q61H",
    interaction_type = "comutation",
    name_suffix = "_TP53",
    genes = c("TP53")
  ),
  tibble(
    cancer = "PAAD",
    allele = "Q61H",
    interaction_type = "exclusivitycomutation",
    name_suffix = "_oncogenes",
    genes = c("TP53", "SMAD4", "ARID1A")
  ),
  tibble(
    cancer = "LUAD",
    allele = "G12C",
    interaction_type = "exclusivity",
    name_suffix = "_FIG5_ABCC9",
    genes = c("ABCC9"),
    keep_gg_onco_proto = TRUE
  ),
  tibble(
    cancer = "LUAD",
    allele = "G12C",
    interaction_type = "exclusivity",
    name_suffix = "_FIG5_CPSF1",
    genes = c("CPSF1"),
    keep_gg_onco_proto = TRUE
  ),
  tibble(
    cancer = "LUAD",
    allele = "G12C",
    interaction_type = "exclusivity",
    name_suffix = "_FIG5_WDFY3",
    genes = c("WDFY3"),
    keep_gg_onco_proto = TRUE
  )
  # tibble(
  #     cancer = "",
  #     allele = "",
  #     interaction_type = "exclusivity comutation",
  #     name_suffix = "",
  #     genes = c()
  # )
) %>%
  group_by(
    cancer, allele, interaction_type, name_suffix, keep_gg_onco_proto
  ) %>%
  summarise(genes = list(genes)) %>%
  ungroup() %>%
  mutate(
    kras_allele = paste0("KRAS_", allele),
    keep_gg_onco_proto = ifelse(is.na(keep_gg_onco_proto),
      FALSE,
      keep_gg_onco_proto
    )
  )

SELECT_DIR <- glue("{GRAPHS_DIR}-select")
reset_graph_directory(SELECT_DIR)

specific_oncoplot_info_tib %>%
  mutate(
    save_name = paste0(
      cancer, "_",
      allele, "_",
      interaction_type,
      "_oncostrip_select",
      name_suffix,
      ".svg"
    )
  ) %>%
  pwalk(
    allele_specificgenes_oncoplot,
    top_n_genes = 15,
    save_dir = SELECT_DIR,
    annotate_kras_allele = TRUE
  )

# Without KRAS allele annotation bar at the bottom.
specific_oncoplot_info_tib %>%
  mutate(
    save_name = paste0(
      cancer, "_",
      allele, "_",
      interaction_type,
      "_oncostrip_select_NO-KRAS-ANNO",
      name_suffix,
      ".svg"
    )
  ) %>%
  pwalk(
    allele_specificgenes_oncoplot,
    top_n_genes = 15,
    save_dir = SELECT_DIR,
    annotate_kras_allele = FALSE
  )


#### ---- Oncoplots for known interactors ---- ####

known_oncoplot_info_tib <- bind_rows(
  tibble(
    cancer = "COAD",
    interaction_type = "comutation",
    name_suffix = "_KNOWN",
    allele = "ALL",
    genes = c("PIK3CA")
  ),
  tibble(
    cancer = "COAD",
    interaction_type = "exclusivity",
    name_suffix = "_KNOWN",
    allele = "ALL",
    genes = c("BRAF", "NRAS")
  ),
  tibble(
    cancer = "LUAD",
    interaction_type = "exclusivity",
    name_suffix = "_KNOWN",
    allele = "ALL",
    genes = c("BRAF", "EGFR")
  ),
  tibble(
    cancer = "LUAD",
    interaction_type = "comutation",
    name_suffix = "_KNOWN",
    allele = "ALL",
    genes = c("STK11")
  ),
  tibble(
    cancer = "MM",
    interaction_type = "exclusivity",
    name_suffix = "_KNOWN",
    allele = "ALL",
    genes = c("NRAS", "BRAF")
  ),
  tibble(
    cancer = "PAAD",
    interaction_type = "exclusivitycomutation",
    name_suffix = "_KNOWN",
    allele = "ALL",
    genes = c("BRAF", "TP53")
  )
) %>%
  group_by(cancer, allele, interaction_type, name_suffix) %>%
  summarise(genes = list(genes)) %>%
  ungroup() %>%
  mutate(kras_allele = paste0("KRAS_", allele))

known_oncoplot_info_tib %>%
  mutate(
    save_name = paste0(
      cancer, "_",
      allele, "_",
      interaction_type,
      "_oncostrip_select",
      name_suffix,
      ".svg"
    )
  ) %>%
  pwalk(
    allele_specificgenes_oncoplot,
    top_n_genes = 15,
    save_dir = SELECT_DIR,
    replace_kras_with_allele = FALSE,
    annotate_kras_allele = TRUE
  )


#### ---- Oncoplots for a priori lists ---- ####
# Oncoplots to accompany the subsetted genetic networks.

gene_sets <- list(
  kegg = unique(kegg_geneset_df$hugo_symbol),
  cgc = unique(filter(cosmic_cgc_df, tier == 1)$hugo_symbol),
  bioid = unique(kras_interactors_bioid_df$hugo_symbol)
)

APRIORI_DIR <- glue("{GRAPHS_DIR}-apriori")
reset_graph_directory(APRIORI_DIR)

for (i in 1:length(gene_sets)) {
  gene_set_name <- names(gene_sets)[[i]]
  gene_set <- unname(unlist(gene_sets[i]))


  cancer_counts_df %>%
    mutate(
      kras_allele = allele,
      allele = str_remove_all(allele, "KRAS_"),
      save_name = paste0(
        cancer, "_",
        allele, "_",
        !!gene_set_name, "_",
        "oncoplot.svg"
      )
    ) %>%
    pwalk(allele_specificgenes_oncoplot,
      genes = gene_set,
      top_n_genes = 10,
      save_dir = APRIORI_DIR,
      annotate_kras_allele = TRUE,
      print_missing_genes = FALSE
    )
}



#### ---- Oncoplots for enriched functions ---- ####

ENRICHED_DIR <- glue("{GRAPHS_DIR}-enriched")
reset_graph_directory(ENRICHED_DIR)

allele_enriched_functions_oncoplot <- function(cancer,
                                               datasource,
                                               term,
                                               mod_term,
                                               allele,
                                               overlap_genes,
                                               ...) {
  allele_specificgenes_oncoplot(
    cancer = cancer,
    kras_allele = paste0("KRAS_", allele),
    genes = overlap_genes,
    interaction_type = "exclusivitycomutation",
    save_name = glue("{cancer}_{allele}_{datasource}_{term}_oncoplot.svg"),
    save_dir = ENRICHED_DIR,
    ignore_genes = NULL,
    top_n_genes = 50,
    replace_kras_with_allele = TRUE,
    annotate_kras_allele = TRUE,
    print_missing_genes = TRUE,
    keep_gg_onco_proto = FALSE,
    ...
  )
}


enrichr_tib %>%
  select(-gene_list) %>%
  unnest(cols = enrichr_res) %>%
  mutate(n_overlap = get_enrichr_overlap_int(overlap)) %>%
  filter(adjusted_p_value < 0.05 & n_overlap >= 3) %>%
  mutate(overlap_genes = str_split(genes, ";")) %>%
  filter(!str_detect(term, !!uninteresting_enrichr_regex)) %>%
  mutate(
    mod_term = standardize_enricher_terms(term),
    mod_term = str_wrap(mod_term, 40),
    mod_term = map2_chr(mod_term, datasource, mod_term_for_datasource)
  ) %>%
  select(-genes) %>%
  filter(n_overlap >= 5) %>%
  pwalk(allele_enriched_functions_oncoplot)
