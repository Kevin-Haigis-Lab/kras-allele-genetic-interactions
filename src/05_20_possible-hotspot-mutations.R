# The possible amino acids from SBS at the hotspots on KRAS.

TABLES_DIR <- "05_20_possible-hotspot-mutations"
reset_table_directory(TABLES_DIR)

save_aa_table <- function(tbl, pos) {
  name <- glue("possible_amino_acids_codon{pos}.tsv")
  write_tsv(tbl, table_path(TABLES_DIR, name))
}


print_possible_mutations <- function(position, codon, ...) {
  cat(paste0(rep("=", 50), collapse = ""), "\n")
  cat(glue("CODON {position} ({codon})"), "\n")
  mut_df <- enumerate_sbs_mutations(codon)
  save_aa_table(mut_df, position)
  all_possible_aa <- sort(unique(mut_df$amino_acid))
  all_possible_aa <- paste(all_possible_aa, collapse = ", ")
  print(mut_df)
  cat("possible amino aicds:", all_possible_aa, "\n")
  cat(paste0(rep("=", 50), collapse = ""), "\n")
  return(all_possible_aa)
}

tibble(
  position = c(12, 13, 61, 146),
  codon = c("GGT", "GGC", "CAA", "GCA"),
  mutants = map2_chr(position, codon, print_possible_mutations)
)
# > # A tibble: 4 x 3
# >   position codon mutants
# >      <dbl> <chr> <chr>
# > 1       12 GGT   A, C, D, G, R, S, V
# > 2       13 GGC   A, C, D, G, R, S, V
# > 3       61 CAA   *, E, H, K, L, P, Q, R
# > 4      146 GCA   A, E, G, P, S, T, V
