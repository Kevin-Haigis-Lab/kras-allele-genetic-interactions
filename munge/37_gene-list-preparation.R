
## Preparation of the gene lists in "data/gene-lists/".


ProjectTemplate::cache("cosmic_cgc_df",
{
    info(logger, "Beginning preparation of COSMIC CGC genes.")
    cosmic_cgc_df <- file.path(
            "data", "gene-lists", "COSMIC-CGC-all_2019-10-10.tsv"
        ) %>%
        read_tsv(col_types = cols()) %>%
        janitor::clean_names() %>%
        rename(hugo_symbol = "gene_symbol") %>%
        mutate(tumour_types_somatic = str_to_lower(tumour_types_somatic)) %>%
        separate_rows(tumour_types_somatic, sep = ",") %>%
        mutate(tumour_types_somatic = str_trim(tumour_types_somatic))

    cosmic_cancers <- cosmic_cgc_df %>%
        filter(!is.na(tumour_types_somatic)) %>%
        u_pull(tumour_types_somatic)

    tumor_type_to_cancer <- list(
        COAD = cosmic_cancers[str_detect(cosmic_cancers, "colon|colorec|^crc$")],
        LUAD = cosmic_cancers[str_detect(cosmic_cancers, "lung|^nsclc$|^luad$")],
        PAAD = cosmic_cancers[str_detect(cosmic_cancers, "^pancre|^pdac$|^paad$")],
        MM = cosmic_cancers[str_detect(cosmic_cancers, "myeloma$|^mm$")]
    ) %>%
        enframe("cancer", "tumour_types_somatic") %>%
        unnest(tumour_types_somatic)

    cosmic_cgc_df <- cosmic_cgc_df %>%
        inner_join(tumor_type_to_cancer, by = "tumour_types_somatic") %>%
        distinct()

    log_rows(logger, cosmic_cgc_df, "cosmic_cgc_df")
    info(logger, "Caching COSMIC CGC data frame.")
    return(cosmic_cgc_df)
})



ProjectTemplate::cache("kegg_geneset_df",
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
            "data", "gsea", "genesets", "c2.cp.kegg.v7.0.symbols.gmt"
        ) %>%
        readgmt::read_gmt(tidy = TRUE) %>%
        filter(str_detect(gene_set, paste0(genesets_to_keep, collapse = "|"))) %>%
        mutate(
            gene_set = str_remove(gene_set, "KEGG_"),
            gene_set = str_replace_us(gene_set),
            gene_set = str_to_sentence(gene_set)
        ) %>%
        dplyr::rename(hugo_symbol = "gene")

    manual_additions <- tribble(
        ~ gene_set, ~ hugo_symbol,
        "KEGG_WNT_SIGNALING_PATHWAY", "AMER1"
    )

    kegg_geneset_df <- bind_rows(kegg_geneset_df, manual_additions)

    log_rows(logger, kegg_geneset_df, "kegg_geneset_df")
    info(logger, "Caching KEGG gene set data frame.")
    return(kegg_geneset_df)
})



ProjectTemplate::cache("kras_interactors_bioid_df",
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



ProjectTemplate::cache("genes_of_interest_df",
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



ProjectTemplate::cache("kea_geneset_df",
{
    kea_path <- file.path("data", "gene-lists", "KEA_2015.txt")
    kea_geneset_df <- readgmt::read_gmt(kea_path, tidy = TRUE) %>%
        filter(gene != "")
    log_rows(logger, kea_geneset_df, "kea_geneset_df")
    info(logger, "Caching data frame of KEA gene set.")
    return(kea_geneset_df)
})



ProjectTemplate::cache("ppiHub_geneset_df",
{
    pppiHub_path <- file.path("data", "gene-lists", "PPI_Hub_Proteins.txt")
    ppiHub_geneset_df <- readgmt::read_gmt(pppiHub_path, tidy = TRUE) %>%
        filter(gene != "")
    log_rows(logger, ppiHub_geneset_df, "ppiHub_geneset_df")
    info(logger, "Caching data frame of PPI Hub gene set.")
    return(ppiHub_geneset_df)
})



ProjectTemplate::cache("chea_geneset_df",
{
    pchea_path <- file.path("data", "gene-lists", "ChEA_2016.txt")
    chea_geneset_df <- readgmt::read_gmt(pchea_path, tidy = TRUE) %>%
        filter(gene != "")
    log_rows(logger, chea_geneset_df, "chea_geneset_df")
    info(logger, "Caching data frame of ChEA gene set.")
    return(chea_geneset_df)
})



ProjectTemplate::cache("tf2dna_tfs",
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



ProjectTemplate::cache("encode_tf_bindingsites",
{
    encode_tf_bindingsites_dir <- file.path(
        "data", "gene-lists", "encode_transcription_factor_targets.gmt"
    )
    encode_tf_bindingsites <- readgmt::read_gmt(encode_tf_bindingsites_dir,
                                                tidy = TRUE)
    log_rows(logger, encode_tf_bindingsites, "encode_tf_bindingsites")
    info(logger, "Caching data frame of ENCODE TF binding sites.")
    return(encode_tf_bindingsites)
})



ProjectTemplate::cache("msigdb_c2_df",
{
    c2_all_path <- file.path("data", "gsea", "genesets", "c2.all.v7.0.symbols.gmt")
    msigdb_c2_df <- readgmt::read_gmt(c2_all_path, tidy = TRUE)
    log_rows(logger, msigdb_c2_df, "msigdb_c2_df")
    info(logger, "Caching data frame of MSigDB C2 gene set.")
    return(msigdb_c2_df)
})



ProjectTemplate::cache("msigdb_hallmark_df",
{
    h_all_path <- file.path("data", "gsea", "genesets", "h.all.v7.0.symbols.gmt")
    msigdb_hallmark_df <- readgmt::read_gmt(h_all_path, tidy = TRUE)
    log_rows(logger, msigdb_hallmark_df, "msigdb_hallmark_df")
    info(logger, "Caching data frame of MSigDB Hallmark gene set.")
    return(msigdb_hallmark_df)
})



ProjectTemplate::cache("mutation_3d_hotspots",
{
    mutation_3d_hotspots_path <- file.path("data", "gene-lists",
                                           "gao_2017_3d-hotspots.xls")

    parse_3d_hotspot_mutations <- function(x) {
        str_split(x, ";") %>%
            unlist() %>%
            str_remove_all("\\([:digit:]+\\)")
    }

    mutation_3d_hotspots <- readxl::read_excel(mutation_3d_hotspots_path,
                                               sheet = "Table S1") %>%
        janitor::clean_names() %>%
        select(-data_available_under_odc_open_database_license_o_db_l,
               -cluster) %>%
        dplyr::rename(hugo_symbol = "gene") %>%
        mutate(
            protein_change = purrr::map(residues_number_mutations,
                                        parse_3d_hotspot_mutations),
            p_value = as.numeric(str_remove_all(p_value, "<"))
        ) %>%
        unnest(protein_change) %>%
        select(hugo_symbol, protein_change, p_value)

    stopifnot(!any(is.na(mutation_3d_hotspots)))

    return(mutation_3d_hotspots)
})
