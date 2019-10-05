
#' Get Hugo identifiers from the DepMap combined name format
#'
#' Extracts the Hugo gene identifiers from a vector of values with the format:
#' "Hugo (Entrez)".
#'
get_hugo_from_depmap_ids <- function(x, split_pattern = " ") {
    y <- unlist(stringr::str_split_fixed(x, split_pattern, n = 2)[, 1])
    return(y)
}


#' Get Entrez identifiers from the DepMap combined name format
#'
#' Extracts the Hugo gene identifiers from a vector of values with the format:
#' "Hugo (Entrez)".
#'
get_entrez_from_depmap_ids <- function(x, convert_to_num = FALSE) {
    y <- stringr::str_extract(x, "(?<=\\()[:digit:]+(?=\\))")
    if (convert_to_num) { y <- as.numeric(y) }
    return(y)
}



# helpful data for identifying KRAS mutations
kras_hotspot_codons <- list(num = c(12, 13, 61, 117, 146))
kras_hotspot_codons$char <- as.character(kras_hotspot_codons$num)
kras_hotspot_codons$regex <- paste0(kras_hotspot_codons$char, collapse = "|")



# coding mutation variant classification terms
coding_mut_var_classes <- c(
    "De_novo_Start_OutOfFrame", "Frame_Shift_Del", "Frame_Shift_Ins",
    "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation",
    "Nonstop_Mutation", "Splice_Site", "Start_Codon_Del", "Start_Codon_Ins",
    "Start_Codon_SNP", "Stop_Codon_Del", "Stop_Codon_Ins"
)



# colors for KRAS alleles
allele_palette <- c(
    "KRAS G12C" = "#0a4f4e",
    "KRAS G12D" = "#427ff5",
    "KRAS G12R" = "#af5fe4",
    "KRAS G12V" = "#f8d147",
    "KRAS Other" = "grey76",
    "KRAS Q61H" = "#fd2c3b",
    "KRAS Q61R" = "#f33bea",
    "KRAS G12A" = "#fb7810",
    "KRAS G12S" = "#f4cacb",
    "KRAS G13C" = "#afe642",
    "KRAS G13D" = "#1c4585",
    "KRAS A146T" = "#859947",
    "KRAS A146V" = "#9f75a7",
    "KRAS Q61L" = "#56eead"
)
allele_palette <- c(allele_palette, "WT" = "grey50")
short_allele_pal <- allele_palette
names(short_allele_pal) <- str_remove_all(names(short_allele_pal), "KRAS ")


# colors for cancers
cancer_palette <- c(
    COAD = "darkolivegreen3",
    LUAD = "black",
    MM = "cornflowerblue",
    PAAD = "brown",
    SKCM = "darksalmon"
)