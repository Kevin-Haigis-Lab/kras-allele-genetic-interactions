# Prepare KEGG pathways into tidygraphs

ProjectTemplate::cache("kegg_pathway_grs", {
  kegg_pathway_data_dir <- file.path("data", "kegg-pathways")
  kegg_pathway_paths <- list.files(kegg_pathway_data_dir,
    pattern = "^hsa",
    full.names = TRUE
  )

  kegg_pathway_grs <- rep(NA, length(kegg_pathway_paths))
  names(kegg_pathway_grs) <- basename(kegg_pathway_paths) %>%
    file_sans_ext() %>%
    str_remove("^hsa[:digit:]+_") %>%
    janitor::make_clean_names()

  for (i in 1:length(kegg_pathway_grs)) {
    path <- kegg_pathway_paths[[i]]
    kegg_pathway_grs[i] <- list(parse_and_annotate_kegg_kgml(path))
  }
  return(kegg_pathway_grs)
})
