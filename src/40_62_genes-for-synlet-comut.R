
pastel_red <- "#FF8FAC"
pastel_blue <- "#4096B3"


#### ---- MASKING ---- ####


mask_pal <- c(short_allele_pal["G12D"],
  TP53 = "#FF9945",
  SMAD4 = "#AACC7C"
)

# All of the groups to make masking group plots for.
masking_hits <- tribble(
  ~cancer, ~allele, ~hugo_symbol, ~other_gene,
  "COAD", "G12D", "STARD9", "TP53",
  "PAAD", "G12D", "EEF1E1", "SMAD4",
  "PAAD", "G12D", "MYBL2", "SMAD4",
  "PAAD", "G12D", "ABI1", "SMAD4",
)
