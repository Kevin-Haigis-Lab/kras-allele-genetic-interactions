# A function that describes every possible amino acid that can be made by
# single base substitutions to a codon.

# The possible DNA nucleotides
DNA_NUCLEOTIDES <- c("A", "C", "G", "T")

# Replace a single nucleotide in x at a position (numeric) with a
# nucleotide (character)
replace_nucleotide <- function(x, at, with = "") {
  Biostrings::replaceAt(x, IRanges::IRanges(at, at), with)
}


# Translate a codon and return the character (single-letter code).
translate_codon <- function(codon) {
  as.character(Biostrings::translate(codon))
}


# Make a table of the possible mutations and amino acid outcomes.
# make_mutation_tibble <- function(muts) {
#     tibble(
#         codon = purrr::map_chr(muts, as.character),
#         amino_acid = purrr::map_chr(muts, translate_codon)
#     )
# }

make_mutation_tibble <- function(new_nuc, pos, codon) {
  aa <- translate_codon(codon)
  tibble(
    nuc = new_nuc,
    pos = pos,
    codon = as.character(codon),
    amino_acid = aa
  )
}


# Enumerate all of the single-base-substitutions possible on a codon.
enumerate_sbs_mutations <- function(codon, silent = FALSE) {
  stopifnot(str_length(codon) == 3)

  codon_bio <- Biostrings::DNAString(codon)
  mut_df <- make_mutation_tibble("-", 0, codon_bio)

  for (N in DNA_NUCLEOTIDES) {
    for (i in seq(1, length(codon_bio))) {
      mut <- replace_nucleotide(codon_bio, i, N)
      df <- make_mutation_tibble(N, i, mut)
      mut_df <- bind_rows(mut_df, df)
    }
  }

  if (!silent) {
    cat(
      glue("=> Num. possible amino acids: {n_distinct(mut_df$amino_acid)}"),
      "\n"
    )
  }
  return(mut_df)
}

# EXAMPLE: codon 12 of KRAS
# enumerate_sbs_mutations("GGT")
# > => Num. possible amino acids: 7
# > # A tibble: 13 x 4
# >    nuc     pos codon amino_acid
# >    <chr> <dbl> <chr> <chr>
# >  1 -         0 GGT   G
# >  2 A         1 AGT   S
# >  3 A         2 GAT   D
# >  4 A         3 GGA   G
# >  5 C         1 CGT   R
# >  6 C         2 GCT   A
# >  7 C         3 GGC   G
# >  8 G         1 GGT   G
# >  9 G         2 GGT   G
# > 10 G         3 GGG   G
# > 11 T         1 TGT   C
# > 12 T         2 GTT   V
# > 13 T         3 GGT   G
