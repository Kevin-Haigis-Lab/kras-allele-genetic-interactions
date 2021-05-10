# Observed vs. predicted KRAS Q61 allele frequency using the trinucleotide
# context.

set.seed(0)

Q61_CODON <- "CAA"
Q61_5NUC_CONTEXT <- ""

GRAPHS_DIR <- "50_17_kras-q61-allele-frequencies"
reset_graph_directory(GRAPHS_DIR)
reset_table_directory(GRAPHS_DIR)


mut_sig_descriptions <- get_mut_sig_descriptions()
artifact_signatures <- get_artifact_signatures()
msg_sigs <- paste0(artifact_signatures, collapse = ", ")
message(glue("The following are artifact signatures: {msg_sigs}"))

# Possible mutations to Q61 codon (ignoring silent and nonsense muts).
possible_q61_alleles <- enumerate_sbs_mutations(Q61_CODON, silent = TRUE) %>%
  filter(amino_acid != "Q") %>%
  filter(amino_acid != "*") %>%
  mutate(kras_allele = as.character(glue("Q61{amino_acid}")))


observed_kras_q61_mutations <- trinucleotide_mutations_df %>%
  filter(hugo_symbol == "KRAS") %>%
  filter(kras_allele != "WT") %>%
  filter(str_detect(kras_allele, "Q61"))

# Assert that all possible mutations were found in the dataset. Can then use the
# mutation's trinucleotide context for further analyses.
stopifnot(
  all(
    possible_q61_alleles$kras_allele %in%
      observed_kras_q61_mutations$amino_acid_change
  )
)

q61_tricontext_mutations <- possible_q61_alleles %>%
  left_join(
    observed_kras_q61_mutations %>%
      select(kras_allele = amino_acid_change, context, tricontext) %>%
      distinct(),
    by = "kras_allele"
  ) %>%
  select(kras_allele, amino_acid, context, tricontext)

knitr::kable(q61_tricontext_mutations)
# |kras_allele |amino_acid |context |tricontext |
# |:-----------|:----------|:-------|:----------|
# |Q61K        |K          |TCA     |T[C>A]A    |
# |Q61P        |P          |TTG     |T[T>G]G    |
# |Q61H        |H          |CTT     |C[T>A]T    |
# |Q61H        |H          |CTT     |C[T>G]T    |
# |Q61E        |E          |TCA     |T[C>G]A    |
# |Q61R        |R          |TTG     |T[T>C]G    |
# |Q61L        |L          |TTG     |T[T>A]G    |
# |Q61H        |H          |CTT     |C[T>A]T    |
# |Q61H        |H          |CTT     |C[T>G]T    |


# Real frequency of Q61 alleles in total populations, just of KRAS mutations,
# and just of KRAS Q61 mutations.

# TODO:


#### ---- Samples with too few mutations ---- ####

MIN_NUM_MUTATION_PER_SAMPLE <- 20

REMOVE_TSB <- trinucleotide_mutations_df %>%
  filter(cancer != "SKCM") %>%
  count(tumor_sample_barcode) %>%
  filter(n < MIN_NUM_MUTATION_PER_SAMPLE) %>%
  pull(tumor_sample_barcode)

message(glue(
  "Removing {length(REMOVE_TSB)} samples because they have too few mutations"
))
