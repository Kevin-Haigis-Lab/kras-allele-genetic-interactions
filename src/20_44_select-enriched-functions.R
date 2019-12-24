# A tibble of the enriched functions to plot for final images.
# These are likely to be used for the figures.

selected_enrichments <- tibble::tribble(
    ~cancer, ~datasource, ~term,
    "COAD", "WikiPathways_2019_Human", "Focal Adhesion",
    "COAD", "Panther_2016", "Hedgehog signaling pathway",
    "COAD", "Reactome_2016", "Beta-catenin phosphorylation cascade",
    "COAD", "Reactome_2016", "Developmental Biology",
    "COAD", "Transcription_Factor_PPIs", "SMAD2",
    "COAD", "KEA_2015", "ATM",
    "COAD", 'PPI_Hub_Proteins', "YWHAZ",
    "COAD", "KEGG_2019_Human", "mTOR signaling pathway",
    "COAD", "KEGG_2019_Human", "Hippo signaling pathway",
    "COAD", "KEGG_2019_Human", "Cellular senescence",
    "COAD", "KEGG_2019_Human", "Apoptosis",
    "COAD", "KEGG_2019_Human", "PI3K-Akt signaling pathway",
    "COAD", "KEGG_2019_Human", "Wnt signaling pathway",
    "LUAD", "WikiPathways_2019_Human", "Focal Adhesion",
    "LUAD", "Panther_2016", "Wnt signaling pathway",
    "LUAD", "GO_Biological_Process_2018", "positive regulation of MAPK cascade",
    "LUAD", "Transcription_Factor_PPIs", "MYC",
    "LUAD", "WikiPathways_2019_Human", "RAC1/PAK1/p38/MMP2 Pathway",
    "LUAD", "WikiPathways_2019_Human", "Pathways Regulating Hippo Signaling",
    "LUAD", "PPI_Hub_Proteins", "PRKACA",
    "LUAD", "PPI_Hub_Proteins", "PRKCA",
    "LUAD", "PPI_Hub_Proteins", "CTNNB1",
    "LUAD", "KEGG_2019_Human", "PI3K-Akt signaling pathway",
    "LUAD", "KEGG_2019_Human", "Apelin signaling pathway",
    "LUAD", "BioCarta_2016", "Chromatin Remodeling by hSWI/SNF ATP-dependent Complexes",
)

cat("Loaded `selected_enrichments`.\n")