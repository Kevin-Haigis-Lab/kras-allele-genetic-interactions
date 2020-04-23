# A tibble of the mutational signatures detected linked to their descriptors.
signature_description_df <- tibble::tribble(
    ~signature, ~description,
           "1",    "1",
           "2",    "2",
           "3",    "3",
           "4",    "4",
           "5",    "5",
           "6",    "MSI",
          "7A",    "UV",
          "7B",    "UV",
          "7C",    "UV",
          "7D",    "UV",
           "8",    "8",
           "9",    "9",
         "10A",    "10",
         "10B",    "10",
         "10C",    "10",
          "13",    "13",
          "14",    "MSI",
          "15",    "MSI",
          "17",    "17",
        "17V2",    "17",
          "18",    "18",
          "20",    "MSI",
          "21",    "MSI",
          "26",    "MSI",
          "30",    "30",
          "38",    "UV",
          "45",    "artifact",
          "55",    "artifact",
         "AR1",    "artifact",
         "AR2",    "artifact",
        "N3V2",    "N3V2"
)


# Palette for all of the mutational signatures.
# Generated with the 'randomcoloR' package.
mutsig_pal <- c(
    "1" = "#F3635E",
    "2" = "#71E989",
    "3" = "#E58E53",
    "4" = "#E480DC",
    "5" = "#56B7E0",
    "6" = "#5BE4DF",
    "7A" = "#92DC87",
    "7B" = "#E1B344",
    "7C" = "#E046D3",
    "7D" = "#E3E39B",
    "8" = "#AAC476",
    "9" = "#DBE93E",
    "10A" = "#5D64DB",
    "10B" = "#7E5B9E",
    "10C" = "#E2A8D9",
    "13" = "#55E4A7",
    "14" = "#C9DE6A",
    "15" = "#CFE1BF",
    "17" = "#50A789",
    "17V2" = "#5737B8",
    "18" = "#E2D1DC",
    "20" = "#72E756",
    "21" = "#E35194",
    "26" = "#B270E4",
    "30" = "#E2B798",
    "38" = "#A56673",
    "45" = "#8B9AE6",
    "55" = "#DA564D",
    "AR1" = "#B3BAE0",
    "AR2" = "#953CE6",
    "N3V2" = "#9A9C59"
)
ggsave_wrapper(
    show_palette(mutsig_pal, font_family = "Arial"),
    plot_path("00_miscellaneous", "mutational-signatures_pal.svg"),
    width = 3.5, height = 2
)

# Palette for all of the mutational signatures.
# Generated with the 'randomcoloR' package.
mutsig_descrpt_pal <- c(
    "1" = "#F3635E",  # VIBRANT
    "2" = "#BA9CBA",  # MUTED
    "3" = "#93CBED",
    "4" = "#86C5B1",  # VIBRANT
    "5" = "#FFAE00",  # VIBRANT
    "6" = "#DBE93E",
    "8" = "#C3E890",  # MUTED
    "9" = "#BF53E0",  # VIBRANT
    "10" = "#C286B5",  # MUTED
    "13" = "#C79271",  # MUTED
    "17" = "#E880C1",  # VIBRANT
    "18" = "#0080FF",  # VIBRANT
    "26" = "#EFA61C",
    "30" = "#E2B798",
    "N3V2" = "#9A9C59",
    "UV" = "#887BC2",  # MUTED
    "MSI" = "#C3768A",  # MUTED
    "artifact" = "gray50"
)
ggsave_wrapper(
    show_palette(mutsig_descrpt_pal, font_family = "Arial"),
    plot_path("00_miscellaneous", "mutational-signature-descriptions_pal.svg"),
    width = 3.5, height = 2
)

# Mutational signature contexts
mutsigs_contexts <- colnames(deconstructSigs::signatures.cosmic)


# Colors from 'MutationalPatterns' R package `COLORS6` object.
# https://rdrr.io/bioc/MutationalPatterns/src/R/MutationalPatterns.R
mutsig_context_group_pal <- c(
    "C>A" = "#2EBAED",
    "C>G" = "#000000",
    "C>T" = "#DE1C14",
    "T>A" = "#D4D2D2",
    "T>C" = "#ADCC54",
    "T>G" = "#F0D0CE"
)

mutsig_context_group_pal <- c(
    "C>A" = "grey60",
    "C>G" = "grey70",
    "C>T" = "grey80",
    "T>A" = "grey55",
    "T>C" = "grey75",
    "T>G" = "grey65"
)





ggsave_wrapper(
    show_palette(mutsig_context_group_pal, font_family = "Arial"),
    plot_path("00_miscellaneous", "mutsig-context-group_pal.svg"),
    width = 3.5, height = 2
)
