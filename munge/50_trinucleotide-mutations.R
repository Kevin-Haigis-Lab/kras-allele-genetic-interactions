# Prepare the trinucleotide mutation data from GM

gm_root_dir <- file.path(
  "/n", "data1", "hms", "dbmi", "park", "gmelloni",
  "cooc_mutex_driver", "data", "cbioportal"
)


ProjectTemplate::cache("trinucleotide_mutations_df", {
  get_data_path <- function(x) file.path("data", "cancer-data", x)

  hypermutants <- get_data_path("hypermut_Aug2019.rds") %>%
    readRDS() %>%
    unlist() %>%
    unique()

  ras_muts <- ras_mutant_tib %>%
    filter(ras == "KRAS") %>%
    mutate(kras_allele = str_remove(ras_allele, "KRAS_")) %>%
    select(cancer, dataset, tumor_sample_barcode, kras_allele)


  # Extract the context from the `tricontext`.
  extract_context <- function(x) {
    nuc1 <- str_sub(x, 1, 1)
    nuc2 <- str_sub(x, 3, 3)
    nuc3 <- str_sub(x, -1, -1)
    return(paste0(nuc1, nuc2, nuc3))
  }

  trinucleotide_mutations_df <- readRDS(file.path(
    gm_root_dir, "paad_luad_coadread_skcm_mm_tricLong_newAug2019.rds"
  )) %>%
    as_tibble() %>%
    janitor::clean_names() %>%
    dplyr::rename(
      tumor_sample_barcode = sample,
      chromosome = chr,
      hugo_symbol = gene_symbol,
      dataset = genetic_profile_id,
      cancer = tumor_type
    ) %>%
    filter(!(tumor_sample_barcode %in% !!double_kras_mutants)) %>%
    mutate(
      mutation_type = str_to_lower(mutation_type),
      mutation_type_hr = unlist(mapping_mutation_types_to_human_readable[mutation_type]),
      is_hypermutant = tumor_sample_barcode %in% hypermutants,
      cancer = str_to_upper(cancer),
      cancer = ifelse(cancer == "COADREAD", "COAD", cancer),
      context = purrr::map_chr(tricontext, extract_context)
    ) %>%
    left_join(ras_muts,
      by = c("cancer", "dataset", "tumor_sample_barcode")
    ) %>%
    mutate(kras_allele = ifelse(is.na(kras_allele), "WT", kras_allele))

  return(trinucleotide_mutations_df)
})



ProjectTemplate::cache("kras_trinucleotide_contexts",
  depends = "trinucleotide_mutations_df",
  {
    kras_trinucleotide_contexts <- trinucleotide_mutations_df %>%
      filter(hugo_symbol == "KRAS") %>%
      filter(amino_acid_change %in% names(short_allele_pal)) %>%
      dplyr::select(amino_acid_change, amino_position, context, tricontext) %>%
      dplyr::rename(
        kras_allele = "amino_acid_change",
        kras_codon = "amino_position"
      ) %>%
      unique() %>%
      mutate(kras_codon = as.numeric(kras_codon)) %>%
      arrange(kras_codon, kras_allele)

    return(kras_trinucleotide_contexts)
  }
)



ProjectTemplate::cache("tricontext_counts_df", {
  tricontext_counts_df <- full_join(
    tibble(
      context = rownames(deconstructSigs::tri.counts.exome),
      exome_count = unlist(deconstructSigs::tri.counts.exome$x)
    ),
    tibble(
      context = rownames(deconstructSigs::tri.counts.genome),
      genome_count = unlist(deconstructSigs::tri.counts.genome$x)
    ),
    by = "context"
  )

  return(tricontext_counts_df)
})
