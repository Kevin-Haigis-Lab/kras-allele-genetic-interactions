
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

# cancer regex
cancer_regex <- "COAD|LUAD|PAAD|MM|SKCM"

# coding mutation variant classification terms
coding_mut_var_classes <- c(
    "De_novo_Start_OutOfFrame", "Frame_Shift_Del", "Frame_Shift_Ins",
    "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation",
    "Nonstop_Mutation", "Splice_Site", "Start_Codon_Del", "Start_Codon_Ins",
    "Start_Codon_SNP", "Stop_Codon_Del", "Stop_Codon_Ins"
)


mapping_mutation_types_to_human_readable <- list(
    "missense_mutation" = "missense mutation",
    "frame_shift_ins" = "frameshift insertion",
    "frame_shift_del" = "frameshift deletion",
    "3'utr" = "3' UTR",
    "intron" = "intronic",
    "in_frame_del" = "in-frame deletion",
    "5'utr" = "5' UTR",
    "5'flank" = "5' flank",
    "silent" = "silent",
    "3'flank" = "3' flank",
    "splice_site" = "splice-site",
    "nonstop_mutation" = "nonstop mutation",
    "nonsense_mutation" = "nonsense mutation",
    "in_frame_ins" = "in-frame insertion",
    "translation_start_site" = "translation start site",
    "splice_region" = "slice region",
    "rna" = "RNA",
    "igr" = "IGR",
    "intronic" = "intronic",
    "3flank" = "3' flank",
    "5flank" = "5' flank",
    "3utr" = "3' UTR",
    "flank" = "flank",
    "5utr" = "5' UTR",
    "inframe_insertion" = "in-frame insertion",
    "frameshift_insertion" = "frameshift insertion",
    "inframe_deletion" = "inframe deletion",
    "frameshift_deletion" = "frameshift deletion",
    "unknown" = "unknown",
    "utr" = "UTR",
    "targeted_region" = "targeted region",
    "de_novo_start_outofframe" = "de novo start (out-of-frame)",
    "start_codon_snp" = "start codon SNV",
    "de_novo_start_inframe" = "de novo start (in-frame)",
    "missense" = "missense mutation"
)


## make color palette: does NOT run on O2, but does on Mac
# set.seed(12)
# n <- 35
# cols <- randomcoloR::distinctColorPalette(n)
# for (i in 1:n) {
#     if (i == 1) { cat("\t") }

#     if (i %% 6 == 0) { cat("\n\t") }
#     cat("\"", cols[[i]], "\", ", sep = "")

#     if (i == n) { cat("\n") }
# }

mutation_pal <- list(
    "#E1B954", "#66E7AE", "#BD9BE3", "#ABC9E0", "#776BE3",
    "#B8887A", "#AAB6A8", "#689658", "#8B569F", "#D26C3F", "#AFE6AF",
    "#73A2DA", "#549C97", "#E9A783", "#EBD4D1", "#74E63D", "#CFC1E8",
    "#E88ED5", "#E44DB5", "#D5E041", "#62E5D7", "#E33E68", "#DC46E7",
    "#D76C8D", "#8942EB", "#5D86E0", "#65CCEA", "#91778D", "#DED5A2",
    "#A8E9E3", "#E1F0D9", "#CC73DF", "#EAA9C9", "#DAED88", "#7BDF7A"
)
names(mutation_pal) <- mapping_mutation_types_to_human_readable



# colors for KRAS alleles
allele_palette <- c(
    "KRAS G12C" = "#0a4f4e",
    "KRAS G12D" = "#427ff5",
    "KRAS G12R" = "#af5fe4",
    "KRAS G12V" = "#f8d147",
    "KRAS G12A" = "#fb7810",
    "KRAS G12S" = "#f4cacb",
    "KRAS G13C" = "#afe642",
    "KRAS G13D" = "#1c4585",
    "KRAS A146T" = "#859947",
    "KRAS A146V" = "#9f75a7",
    "KRAS Q61L" = "#56eead",
    "KRAS Q61K" = "#E76707",
    "KRAS Q61H" = "#fd2c3b",
    "KRAS Q61R" = "#f33bea",
    "KRAS K117N" = "#E707C8",
    "KRAS Other" = "grey75"
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


# colors for comutation vs mutual exclusivity
comut_mutex_pal <- c(
    comutation = "lightskyblue",
    exclusivity = "palegreen3"
)
