# Cancer specific helpers.

# Specific list of known and prominent oncogenes for each cancer.
cancer_oncogenes <- list(
    COAD = c(
        "APC", "TP53"
    ),
    LUAD = c(
        "EGFR", "TP53", "STK11"
    ),
    MM = c(
        "NRAS", "TP53", "DIS3", "FAM46C", "BRAF", "TRAF3", "PRDM1", "CYLD",
        "RB1", "ACTG1"
    ),
    PAAD = c(
        "TP53", "CDKN2A", "SMAD4"
    )
)


FRACTION_OF_TISSUE_THAT_ARE_CANCER <- list(
    "COAD" = 0.96,  # ACS: https://www.cancer.org/cancer/colon-rectal-cancer/about/what-is-colorectal-cancer.html
    "LUAD" = 0.50,  # PMID: 25822850
      "MM" = 1.00,
    "PAAD" = 0.95,   # ACS: https://www.cancer.org/cancer/pancreatic-cancer/about/what-is-pancreatic-cancer.html
    "SKCM" = 1.00
)


# Regular expressions to filter `cosmic_cgc_df$tumour_types_somatic`.
cgc_cancer_regex <- list(
    COAD = "colo|rect",
    LUAD = "lung|NSCLC",
    MM = "myeloma",
    PAAD = "pancrea"
)

# Get the genes labeled by the CGC as oncogenes or TS for a cancer.
# Set `full_df=TRUE` if the entire COSMIC data frame is desired.
get_cgc_genes <- function(cancer, tier = c(1, 2), full_df = FALSE) {
    if (!cancer %in% names(cgc_cancer_regex)) {
        stop(glue("No cancer '{cancer}'"))
    }

    reg_exp <- cgc_cancer_regex[[cancer]]

    df <- cosmic_cgc_df %>%
        filter(
            str_detect(tumour_types_somatic, !!reg_exp) |
            str_detect(tumour_types_germline, !!reg_exp)
        )

    if (full_df) {
        return(df)
    }
    return(unique(df$hugo_symbol))
}
