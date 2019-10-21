
## Preparation of the gene lists in "data/gene-lists/".


cache("cosmic_cgc_df",
{
    info(logger, "Beggining preparation of COSMIC CGC genes.")
    cosmic_cgc_df <- file.path(
            "data", "gene-lists", "COSMIC-CGC-all_2019-10-10.tsv"
        ) %>%
        read_tsv(col_types = cols()) %>%
        janitor::clean_names() %>%
        dplyr::rename(hugo_symbol = "gene_symbol")

    cosmic_cancers <- cosmic_cgc_df %>%
        pull(tumour_types_somatic) %>%
        str_split(",") %>%
        unlist() %>%
        str_remove_all(" ") %>%
        unique()

    relevant_cosmic_cancers <- cosmic_cancers[str_detect(cosmic_cancers, "colon|lung|NSCLC|pancreas|myeloma")]
    relevant_cosmic_cancers <- relevant_cosmic_cancers[!is.na(relevant_cosmic_cancers)]
    regex_cosmic_cancers <- paste0(relevant_cosmic_cancers, collapse = "|")

    cosmic_cgc_df <- cosmic_cgc_df %>%
        filter(
            str_detect(tumour_types_somatic, regex_cosmic_cancers) |
            str_detect(tumour_types_germline, regex_cosmic_cancers)
        )

    log_rows(logger, cosmic_cgc_df, "cosmic_cgc_df")
    info(logger, "Caching COSMIC CGC data frame.")
    return(cosmic_cgc_df)
})



cache("kegg_geneset_df",
{

    genesets_to_keep <- c(
        "KEGG_APOPTOSIS",
        "KEGG_CELL_CYCLE",
        "KEGG_ERBB_SIGNALING_PATHWAY",
        "KEGG_HEDGEHOG_SIGNALING_PATHWAY",
        "KEGG_MAPK_SIGNALING_PATHWAY",
        "KEGG_MTOR_SIGNALING_PATHWAY",
        "KEGG_NOTCH_SIGNALING_PATHWAY",
        "KEGG_P53_SIGNALING_PATHWAY",
        "KEGG_VEGF_SIGNALING_PATHWAY",
        "KEGG_WNT_SIGNALING_PATHWAY"
    )

    kegg_geneset_df <- file.path(
            "data", "gene-lists", "c2.cp.kegg.v7.0.symbols.gmt"
        ) %>%
        readgmt::read_gmt(tidy = TRUE) %>%
        filter(str_detect(gene_set, paste0(genesets_to_keep, collapse = "|"))) %>%
        mutate(
            gene_set = str_remove(gene_set, "KEGG_"),
            gene_set = str_replace_us(gene_set),
            gene_set = str_to_sentence(gene_set)
        ) %>%
        dplyr::rename(hugo_symbol = "gene")

    log_rows(logger, kegg_geneset_df, "kegg_geneset_df")
    info(logger, "Caching KEGG gene set data frame.")
    return(kegg_geneset_df)
})



cache("kras_interactors_bioid_df",
{
    kras_interactors_bioid_df <- file.path(
            "data", "gene-lists", "Kovalski-et-al-2019_BioID_REFORMATTED.xlsx"
        ) %>%
        readxl::read_xlsx() %>%
        janitor::clean_names() %>%
        gather("key", "value", -gene) %>%
        filter(key == "kras_wt_hit" & value == 1) %>%
        dplyr::rename(hugo_symbol = "gene")

    log_rows(logger, kras_interactors_bioid_df, "kras_interactors_bioid_df")
    info(logger, "Caching data frame of KRAS interactors from BioID study.")
    return(kras_interactors_bioid_df)
})



cache("genes_of_interest_df",
      depends = c("cosmic_cgc_df",
                  "kegg_geneset_df",
                  "kras_interactors_bioid_df"),
{
    genes_of_interest_df <- bind_rows(
        {
            cosmic_cgc_df %>%
                select(hugo_symbol) %>%
                add_column(source = "CGC") %>%
                unique()
        },
        {
            kegg_geneset_df %>%
                dplyr::rename(source = "gene_set") %>%
                unique()
        },
        {
            kras_interactors_bioid_df %>%
                select(hugo_symbol) %>%
                add_column(source = "BioID") %>%
                unique()
        }
    )

    log_rows(logger, genes_of_interest_df, "genes_of_interest_df")
    info(logger, "Caching data frame of list of genes of interest.")
    return(genes_of_interest_df)
})



cache("kea_geneset_df",
{
    kea_path <- file.path("data", "gene-lists", "KEA_2015.txt")
    kea_geneset_df <- readgmt::read_gmt(kea_path, tidy = TRUE) %>%
        filter(gene != "")
    log_rows(logger, kea_geneset_df, "kea_geneset_df")
    info(logger, "Caching data frame of KEA gene set.")
    return(kea_geneset_df)
})



cache("ppiHub_geneset_df",
{
    pppiHub_path <- file.path("data", "gene-lists", "PPI_Hub_Proteins.txt")
    ppiHub_geneset_df <- readgmt::read_gmt(kea_path, tidy = TRUE) %>%
        filter(gene != "")
    log_rows(logger, ppiHub_geneset_df, "ppiHub_geneset_df")
    info(logger, "Caching data frame of PPI Hub gene set.")
    return(ppiHub_geneset_df)
})




cache("chea_geneset_df",
{
    pchea_path <- file.path("data", "gene-lists", "ChEA_2016.txt")
    chea_geneset_df <- readgmt::read_gmt(kea_path, tidy = TRUE) %>%
        filter(gene != "")
    log_rows(logger, chea_geneset_df, "chea_geneset_df")
    info(logger, "Caching data frame of ChEA gene set.")
    return(chea_geneset_df)
})



cache("tf2dna_tfs",
{
    tf2dna_dir <- file.path("data", "gene-lists", "human_tf2dna_matrices_symbols")
    tf2dna_tfs <- list.files(tf2dna_dir) %>%
        unlist() %>%
        basename() %>%
        str_remove_all("\\.mat$") %>%
        unique()
    info(logger, "Caching transcription factor list from TF2DNA.")
    return(tf2dna_tfs)
})




cache("msigdb_c2_df",
{
    c2_all_path <- file.path("data", "gene-lists", "c2.all.v7.0.symbols.gmt")
    msigdb_c2_df <- readgmt::read_gmt(c2_all_path, tidy = TRUE)
    log_rows(logger, msigdb_c2_df, "msigdb_c2_df")
    info(logger, "Caching data frame of MSigDB C2 gene set.")
    return(msigdb_c2_df)
})
