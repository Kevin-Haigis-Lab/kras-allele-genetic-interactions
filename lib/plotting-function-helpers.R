
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
