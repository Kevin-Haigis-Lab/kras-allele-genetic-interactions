# Miscellaneous things for KH.


#### ---- KRAS & NRAS comutation in COAD ---- ####

kras_nras_comuts <- cancer_full_coding_muts_df %>%
  filter(cancer == "COAD") %>%
  filter(hugo_symbol %in% c("KRAS", "NRAS")) %>%
  select(tumor_sample_barcode, hugo_symbol, amino_acid_change) %>%
  pivot_wider(
    id_cols = c(tumor_sample_barcode),
    names_from = hugo_symbol,
    values_from = amino_acid_change,
    values_fn = function(x) {
      x[x == ""] <- "(non-coding)"
      if (is.null(x)) {
        return("WT")
      } else {
        return(paste(x, collapse = ", "))
      }
    }
  ) %>%
  mutate(
    KRAS = ifelse(is.na(KRAS), "WT", KRAS),
    NRAS = ifelse(is.na(NRAS), "WT", NRAS)
  ) %>%
  filter(KRAS != "WT" & NRAS != "WT")

kras_nras_comuts %>%
  knitr::kable(format = "simple")

num_samples <- n_distinct(cancer_full_coding_muts_df$tumor_sample_barcode)
num_kras_muts <- cancer_full_coding_muts_df %>%
  filter(ras_allele != "WT") %>%
  pull(tumor_sample_barcode) %>%
  n_distinct()
num_nras_muts <- cancer_full_coding_muts_df %>%
  filter(hugo_symbol == "NRAS") %>%
  pull(tumor_sample_barcode) %>%
  n_distinct()

cat(
  paste(
    glue("          total number of tumor samples: {num_samples}"),
    glue("           total number of KRAS mutants: {num_kras_muts}"),
    glue("           total number of NRAS mutants: {num_nras_muts}"),
    glue("total number of KRAS / NRAS double muts: {nrow(kras_nras_comuts)}"),
    sep = "\n"
  ),
  "\n"
)

################################################################################
################################################################################
