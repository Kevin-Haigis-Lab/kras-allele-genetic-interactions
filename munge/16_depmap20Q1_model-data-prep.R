# Prepare the data to be used to model depletion effects from the KRAS allele


#### ---- Double KRAS mutants ---- ####

cache("ccle_kras_double_muts", depends = "ccle_kras_muts", {
  ccle_kras_double_muts <- ccle_kras_muts %>%
    filter(allele != "other") %>%
    distinct(dep_map_id, allele) %>%
    add_count(dep_map_id, name = "total_kras_muts") %>%
    filter(total_kras_muts > 1)
})



#### ---- MAPK mutants ---- ####

cache("mapk_muts", depends = "ccle_mutations_dmg", {
  mapk_muts <- ccle_mutations_dmg %>%
    filter(
      (hugo_symbol == "NRAS" &
        str_detect(protein_change, "G12|G13|Q61")) |
        (hugo_symbol == "BRAF" & str_detect(protein_change, "V600")) |
        (hugo_symbol == "EGFR" & is_cosmic_hotspot)
    ) %>%
    distinct()
  return(mapk_muts)
})



#### ---- Gene is altered ---- ####

# A table of gene deleted ("homo_del") in cell lines.
# These data points should be ignored when modeling because the gene is not
# present to be knocked out.
gene_is_deleted <- ccle_copy_number %>%
  mutate(is_deleted = copy_number_label == "homo_del") %>%
  select(dep_map_id, hugo_symbol, is_deleted)


# A table of mutated genes in the cell lines.
gene_is_mutated <- ccle_mutations_dmg %>%
  mutate(is_mutated = TRUE) %>%
  select(dep_map_id, hugo_symbol, is_mutated)



#### ---- Genes that are not expressed in normal or cancer tissue ---- ####

# Genes that are confidently not expressed in tumors.
# These should be removed from their respective cancers.
unexpressed_gene_tibble <- enframe(confidently_unexpressed_genes) %>%
  rename(
    cancer = "name",
    hugo_symbol = "value"
  ) %>%
  mutate(
    cancer = str_to_upper(cancer),
    is_unexpressed = TRUE
  ) %>%
  unnest(hugo_symbol)

# Print out the number of genes to remove
knitr::kable(table(unexpressed_gene_tibble$cancer))



#### ---- Cell line information ---- ####

model_cell_line_info <- ccle_cell_lines %>%
  filter(!is.na(cancer)) %>%
  left_join(ccle_kras_muts, by = "dep_map_id") %>%
  select(dep_map_id, cancer,
    kras_allele = allele,
    kras_codon = codon,
    kras_cn = copy_number
  ) %>%
  distinct() %>%
  mutate(
    kras_allele = ifelse(is.na(kras_allele), "WT", kras_allele),
    kras_codon = ifelse(is.na(kras_codon), "WT", kras_codon),
    kras_cn = ifelse(is.na(kras_cn), 2, kras_cn)
  )



#### ---- Make a data frame to use for modeling using DepMap data ---- ####


cache("depmap_modelling_df",
  depends = c(
    "gene_effect", "ccle_kras_double_muts", "mapk_muts",
    "unexpressed_gene_tibble", "model_cell_line_info",
    "gene_is_deleted", "gene_is_mutated", "ccle_expression"
  ),
  {
    depmap_modelling_df <- gene_effect %>%
      filter(!dep_map_id %in% c(
        !!ccle_kras_double_muts$dep_map_id,
        !!mapk_muts$dep_map_id
      )) %>%
      inner_join(model_cell_line_info, by = "dep_map_id") %>%
      anti_join(unexpressed_gene_tibble, by = c("cancer", "hugo_symbol")) %>%
      left_join(gene_is_deleted, by = c("dep_map_id", "hugo_symbol")) %>%
      left_join(gene_is_mutated, by = c("dep_map_id", "hugo_symbol")) %>%
      mutate(
        is_deleted = ifelse(is.na(is_deleted), FALSE, is_deleted),
        is_mutated = ifelse(is.na(is_mutated), FALSE, is_mutated)
      ) %>%
      left_join(select(ccle_expression, -entrez_id),
        by = c("dep_map_id", "hugo_symbol")
      )

    # Remove genes if any of the cell lines in the cancer have an
    #   `NA` gene effect.
    depmap_modelling_df %<>%
      group_by(cancer, hugo_symbol) %>%
      filter(!any(is.na(gene_effect))) %>%
      ungroup()

    # Make sure each row is unique.
    depmap_modelling_df %<>% distinct()

    x <- depmap_modelling_df %>%
      count(dep_map_id, hugo_symbol) %>%
      filter(n > 1)
    stopifnot(nrow(x) == 0)

    return(depmap_modelling_df)
  }
)



#### ---- Gene effect comparisons of non/essential genes ---- ####

cache("essential_gene_effect_summary",
  depends = c("depmap_modelling_df", "essentiality_tib"),
  {
    essential_gene_effect_summary <- depmap_modelling_df %>%
      inner_join(essentiality_tib %>% select(-entrez_id),
        by = "hugo_symbol"
      ) %>%
      filter(label %in% c("common_essential", "nonessential")) %>%
      mutate(label = str_remove(label, "common_")) %>%
      group_by(cancer, label) %>%
      summarise(
        avg = mean(gene_effect),
        sd = sd(gene_effect)
      ) %>%
      ungroup() %>%
      pivot_wider(cancer,
        names_from = label,
        values_from = c(avg, sd)
      ) %>%
      mutate(cut_1sd_from_nonessential = avg_nonessential - sd_nonessential)
    return(essential_gene_effect_summary)
  }
)
