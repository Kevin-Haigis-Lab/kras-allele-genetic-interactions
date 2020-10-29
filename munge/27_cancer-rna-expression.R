## Prepare the RNA expression data from human tumor samples.

cache("tcga_coad_rna", {
  dir_path <- file.path("data", "cbioportal", "coad_tcga")
  data_path <- file.path(dir_path, "data_RNA_Seq_v2_expression_median.txt")

  ras_muts <- cancer_full_coding_muts_df %>%
    filter(dataset == "coadread_tcga_pan_can_atlas_2018") %>%
    distinct(tumor_sample_barcode, ras_allele, is_hypermutant)

  col_types <- cols(
    .default = col_double(),
    Hugo_Symbol = col_character(),
    Entrez_Gene_Id = col_character()
  )

  tcga_coad_rna <- read_tsv(
    file = data_path,
    col_types = col_types,
    progress = FALSE
  ) %.% {
    pivot_longer(
      -c("Hugo_Symbol", "Entrez_Gene_Id"),
      names_to = "tumor_specimen_barcode",
      values_to = "rna_expr"
    )
    janitor::clean_names()
    mutate(tumor_sample_barcode = str_remove(tumor_specimen_barcode, "-[:digit:]+$"))
    filter()
    inner_join(ras_muts, by = "tumor_sample_barcode")
    select(
      hugo_symbol, tumor_sample_barcode, rna_expr, ras_allele, is_hypermutant
    )
  }
})

