# 1. Export a list of non-allele-specific KRAS comutation interactions in CRC.
# 2. Compare expression levels of DUSP genes between KRAS mutant CRC.


GRAPHS_DIR <- "90_44_yijang_dusp-analysis"
TABLES_DIR <- GRAPHS_DIR

reset_table_directory(TABLES_DIR)

#### ---- Comutation interaction list ---- ####


nonallele_specific_increased_comutation_df %.% {
  filter(cancer == "COAD" & hugo_symbol != "KRAS")
  filter(p_value < 0.01)
  arrange(p_value)
  mutate(
    geneWT_krasWT = map_dbl(comut_ct_tbl, ~ .x[1, 1]),
    geneMut_krasWT = map_dbl(comut_ct_tbl, ~ .x[2, 1]),
    geneWT_krasMut = map_dbl(comut_ct_tbl, ~ .x[1, 2]),
    geneMut_krasMut = map_dbl(comut_ct_tbl, ~ .x[2, 2]),
  )
  select(hugo_symbol, p_value, odds_ratio, tidyselect::starts_with("gene"))
} %>%
  write_tsv(table_path(TABLES_DIR, "comutation-list.tsv"))



#### ---- RNA expression data ---- ####

tcga_coad_rna %>%
  filter(str_detect(hugo_symbol, "DUSP")) %>%
  saveRDS("coad-dusp-rna_expression.rds")