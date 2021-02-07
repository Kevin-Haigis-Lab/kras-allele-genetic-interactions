
# Make row and column indices after arranging by a column.
make_rows_cols_idx <- function(df, arrange_col, num_cols) {
  df %>%
    arrange(-{{ arrange_col }}) %>%
    mutate(row_idx = ceiling(row_number() / num_cols)) %>%
    group_by(row_idx) %>%
    arrange({{ arrange_col }}) %>%
    mutate(col_idx = row_number()) %>%
    ungroup()
}

# Assign p-value stars.
assign_stars <- function(pval) {
  case_when(
    pval < 0.001 ~ "***",
    pval < 0.01 ~ "**",
    pval < 0.05 ~ "*",
    TRUE ~ NA_character_
  )
}

# Modify terms from a linear model to add italics to gene names.
italicize_genes_in_terms <- function(terms, not_gene_terms = c("kras_allele")) {
  italicize_interaction_term <- function(s) {
    str_split_fixed(s, ":", 2) %>%
      apply(1, function(x) {
        paste0("*", x, "*", collapse = ":")
      }) %>%
      unlist()
  }

  case_when(
    terms %in% !!not_gene_terms ~ terms,
    str_detect(terms, ":") ~ italicize_interaction_term(terms),
    TRUE ~ paste0("*", terms, "*")
  )
}
