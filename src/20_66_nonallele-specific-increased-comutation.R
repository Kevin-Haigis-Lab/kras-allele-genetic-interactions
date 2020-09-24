# Run the non-allele-specific comutation between KRAS and all other genes.



target_panel_genes <- cancer_full_coding_muts_df %>%
  filter(!target %in% c("exome", "genome")) %>%
  filter(cancer != "SKCM") %>%
  distinct(cancer, hugo_symbol, target)

ts_kras_alleles <- cancer_full_coding_muts_df %>%
  filter(cancer != "SKCM") %>%
  distinct(cancer, tumor_sample_barcode, ras_allele)


# Returns a vector of panels with the gene on them.
get_panels_with_gene <- function(cancer, hugo_symbol) {
  target_panel_genes %>%
    filter(cancer == !!cancer & hugo_symbol == !!hugo_symbol) %>%
    u_pull(target)
}


# Make a 2x2 table from logical vectors.
make_lgl_table <- function(x, y) {
  f <- function(a) {
    factor(a, levels = c("FALSE", "TRUE"))
  }
  table(f(x), f(y))
}


# Make a 2x2 contingency table of comutation.
make_comutation_table <- function(cancer, hugo_symbol, mut_df) {
  if (!exists("target_panel_genes")) {
    stop("Need `target_panel_genes` data frame to exist.")
  }
  if (!exists("ts_kras_alleles")) {
    stop("Need `ts_kras_alleles` data frame to exist.")
  }
  panels_to_include <- get_panels_with_gene(cancer, hugo_symbol)
  mut_df %>%
    filter(cancer == !!cancer) %>%
    filter(target %in% c("genome", "exome", panels_to_include)) %>%
    group_by(tumor_sample_barcode) %>%
    summarise(
      is_mutant = any(hugo_symbol == !!hugo_symbol),
      is_kras = unique(ras_allele != "WT")
    ) %>%
    ungroup() %$%
    make_lgl_table(is_mutant, is_kras)
}


# Run the one-sided Fisher's test of association on a contingency table.
run_fisher_greater_test <- function(ct_tbl) {
  fisher.test(ct_tbl, alternative = "greater") %>%
    broom::tidy() %>%
    janitor::clean_names() %>%
    select(odds_ratio = estimate, p_value, conf_low, conf_high)
}

library(tictoc)
tic("Non-allele-specific comutation analysis")
cache(
  "nonallele_specific_increased_comutation_df",
  depends = c("cancer_full_coding_muts_df"),
  {
    nonallele_specific_increased_comutation_df <- cancer_full_coding_muts_df %>%
      distinct(cancer, hugo_symbol) %>%
      mutate(
        comut_ct_tbl = map2(
          cancer,
          hugo_symbol,
          make_comutation_table,
          mut_df = cancer_full_coding_muts_df
        ),
        fish_grtr_res = map(comut_ct_tbl, run_fisher_greater_test)
      ) %>%
      unnest(fish_grtr_res)

    return(nonallele_specific_increased_comutation_df)
  }
)
toc()
# > Non-allele-specific comutation analysis: 17654.185 sec elapsed
