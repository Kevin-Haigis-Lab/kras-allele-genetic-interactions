################################
## Gene Expression per Tissue ##
################################
# Creates files of tissue-specific gene expression

# normal and cancer tissue gene expression
    # normal tissue databases: GTEx and The Human Protein Atlas (RNA-seq)
    # cancer tissue databases: TCGA (from cBioPortal)

#### ---- Normal Tissue ---- ####
db_dir <- file.path("data", "tissue-gene-expression")
gtex_data_dir <- file.path(db_dir, "GTEx_gene_tpm.txt")
gtex_info_dir <- file.path(db_dir, "GTEx_v7_sampleattributes.txt")
hpa_data_dir <- file.path(db_dir, "HumanProteinAtlas_tissueRNA.txt")

# read in info on each patient
gtex_info <- read_tsv(gtex_info_dir, col_types = cols()) %>%
    select(SAMPID, SMTS) %>%
    dplyr::rename(sampleid = SAMPID,
                  tissue = SMTS) %>%
    filter(tissue %in% c("Colon", "Lung", "Pancreas", "Skin",
                         "Blood", "Bone Marrow")) %>%
    unique()

# a chunk function to melt GTEx info and join with `gtex_info`
f <- function(x, pos) {
    y1 <- x[, colnames(x) %in% c("Name", "Description", gtex_info$sampleid)]
    y2 <- y1 %>%
        gather(key = "sampleid", value = "tpm", -Name, -Description) %>%
        dplyr::rename(ensembl_gene_id = "Name") %>%
        left_join(gtex_info, by = "sampleid")
    return(y2)
}

# GTEX
cache("gtex_summary_expr",
{
    # read in the GTEx data file in chunks
    gtex_data <- read_tsv_chunked(gtex_data_dir,
                                  DataFrameCallback$new(f),
                                  chunk_size = 1E4,
                                  skip = 2)

    # summarise GTEx data into a single value (median and mean)
    gtex_summary_expr <- gtex_data %>%
        mutate(tissue = str_to_lower(tissue)) %>%
        select(Description, sampleid, tpm, tissue) %>%
        group_by(tissue, Description) %>%
        summarise(mid_expr = median(tpm),
                  avg_expr = mean(tpm),
                  num_samples = n_distinct(sampleid)) %>%
        ungroup()
    colnames(gtex_summary_expr) <- c("tissue", "gene", "mid_expr",
                            "avg_expr", "num_samples")
    gtex_summary_expr %<>%
        select(tissue, gene, mid_expr, num_samples) %>%
        dplyr::rename(GTEx_expr = mid_expr,
                      GTEx_num_samples = num_samples) %>%
        mutate(tissue = str_to_lower(tissue))
    return(gtex_summary_expr)
})

# HPA
cache("hpa_expr",
{
    hpa_expr <- read_tsv(hpa_data_dir, col_types = cols()) %>%
        filter(Sample %in% c("colon", "lung", "pancreas", "rectum", "skin",
                             "bone marrow")) %>%
        select(-c("Unit", "Gene"))
    colnames(hpa_expr) <- c("gene", "tissue", "expr")
    hpa_expr %<>%
        select(tissue, gene, expr) %>%
        dplyr::rename(HPA_expr = expr)
    return(hpa_expr)
})

# one data frame with all of the information
cache("normal_tissue_expr",
    depends = c("hpa_expr", "gtex_summary_expr"),
{
    normal_tissue_expr <- full_join(hpa_expr, gtex_summary_expr,
                                    by = c("tissue", "gene"))
    return(normal_tissue_expr)
})

#### ---- Cancer Tissue (TCGA) ---- ####
# paths to the cancer RNA-seq data
raw_data_dir <- file.path(".", "data", "tcga")
coad_rna_dir <- file.path(
    raw_data_dir,
    "coad_tcga_data_expression_merged.txt"
)
luad_rna_dir <- file.path(
    raw_data_dir,
    "luad_tcga_data_expression_merged.txt"
)
paad_rna_dir <- file.path(
    raw_data_dir,
    "paad_tcga_data_RNA_Seq_v2_expression_median.txt"
)
skcm_rna_dir <- file.path(
    raw_data_dir,
    "skcm_tcga_data_RNA_Seq_v2_expression_median.txt"
)
mm_rna_dir <- file.path(
    raw_data_dir,
    "mm_mmrf_MMRF_CoMMpass_IA12a_E74GTF_HtSeq_Gene_Counts.txt"
)

cache("cancer_rna_tib",
{
    # COAD
    coad_rna <- read_tsv(coad_rna_dir, col_types = cols()) %>%
        filter(!is.na(Hugo_Symbol)) %>%
        gather(
            key = "tumor_sample_barcode", value = "rna_expr",
            -Hugo_Symbol, -Entrez_Gene_Id
        ) %>%
        add_column(cancer = "coad") %>%
        select(-Entrez_Gene_Id) %>%
        janitor::clean_names()

    # LUAD
    luad_rna <- read_tsv(luad_rna_dir, col_types = cols()) %>%
        filter(!is.na(Hugo_Symbol)) %>%
        gather(
            key = "tumor_sample_barcode", value = "rna_expr",
            -Hugo_Symbol, -Entrez_Gene_Id
        ) %>%
        add_column(cancer = "luad") %>%
        select(-Entrez_Gene_Id) %>%
        janitor::clean_names()

    # PAAD
    paad_rna <- read_tsv(paad_rna_dir, col_types = cols()) %>%
        filter(!is.na(Hugo_Symbol)) %>%
        gather(
            key = "tumor_sample_barcode", value = "rna_expr",
            -Hugo_Symbol, -Entrez_Gene_Id
        ) %>%
        add_column(cancer = "paad") %>%
        select(-Entrez_Gene_Id) %>%
        janitor::clean_names()

    # SKCM
    skcm_rna <- read_tsv(skcm_rna_dir, col_types = cols()) %>%
        filter(!is.na(Hugo_Symbol)) %>%
        gather(
            key = "tumor_sample_barcode", value = "rna_expr",
            -Hugo_Symbol, -Entrez_Gene_Id
        ) %>%
        add_column(cancer = "skcm") %>%
        select(-Entrez_Gene_Id) %>%
        janitor::clean_names()

    # MM
    mm_rna <- read_tsv(mm_rna_dir, col_types = cols()) %>%
        filter(!is.na(GENE_ID)) %>%
        gather(
            key = "tumor_sample_barcode", value = "rna_expr",
            -GENE_ID) %>%
        add_column(cancer = "mm") %>%
        janitor::clean_names()


    mm_rna <- left_join(
            mm_rna, gene_name_tib,
            by = c("gene_id" = "ensembl_id")
        ) %>%
        select(symbol, tumor_sample_barcode, rna_expr, cancer) %>%
        filter(!is.na(symbol)) %>%
        dplyr::rename(hugo_symbol = "symbol")

    # merge all of the cancers together
    cancer_rna_tib <- bind_rows(
        coad_rna,
        luad_rna,
        paad_rna,
        skcm_rna,
        mm_rna
    )

    # summarise the cancer RNA-seq values (median and mean)
    cancer_rna_tib %<>%
        mutate(rna_expr = ifelse(is.na(rna_expr), 0, rna_expr)) %>%
        group_by(cancer, hugo_symbol) %>%
        summarise(mid_expr = median(rna_expr),
                  avg_expr = mean(rna_expr),
                  num_samples = n_distinct(tumor_sample_barcode)) %>%
        ungroup()

    return(cancer_rna_tib)
})



#### ---- Tissue Expression Gene Lists ---- ####

# return the genes expressed AT OR ABOVE the `min_expr` level
get_tissue_genes <- function(cancer_name, tissue_name, min_expr = 1) {
    tiss <- normal_tissue_expr %>%
        filter(tissue %in% !!tissue_name) %>%
        filter(HPA_expr >= !!min_expr | GTEx_expr >= !!min_expr) %>%
        pull(gene) %>%
        unique()
    canc <- cancer_rna_tib %>%
        filter(cancer %in% !!cancer_name) %>%
        filter(mid_expr >= !!min_expr) %>%
        pull(hugo_symbol) %>%
        unique()
    genes <- unique(c(tiss, canc))
    return(genes)
}

# return the genes expressed BELOW the `max_expr` level
get_tissue_genes_with_max_expr <- function(cancer_name, tissue_name, max_expr = 1) {
    tiss <- normal_tissue_expr %>%
        filter(tissue %in% !!tissue_name) %>%
        filter(HPA_expr < !!max_expr | GTEx_expr < !!max_expr) %>%
        pull(gene) %>%
        unique()
    canc <- cancer_rna_tib %>%
        filter(cancer %in% !!cancer_name) %>%
        filter(mid_expr < !!max_expr) %>%
        pull(hugo_symbol) %>%
        unique()
    genes <- unique(c(tiss, canc))
    return(genes)
}


cache("confidently_unexpressed_genes", {
    confidently_unexpressed_genes <- list(
        coad = get_tissue_genes_with_max_expr("coad", c("rectum", "colon"), max_expr = 1),
        luad = get_tissue_genes_with_max_expr("luad", "lung", max_expr = 1),
        paad = get_tissue_genes_with_max_expr("paad", "pancreas", max_expr = 1),
        skcm = get_tissue_genes_with_max_expr("skcm", "skin", max_expr = 1),
        mm = get_tissue_genes_with_max_expr("mm", c("blood", "bone marrow"), max_expr = 1)
    )
})
