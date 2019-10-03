################################
## Gene Expression per Tissue ##
################################
# Creates files of tissue-specific gene expression


# library(biomaRt)

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

stop("JHC: no done with munging script")

#### ---- Tissue Expression Gene Lists ---- ####

# normal_expr <- readRDS("normal_tissue_expr.tib")
# cancer_expr <- readRDS("cancer_tissue_expr.tib")

# return the genes expressed at the `min_expr` level
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

# make lists of expressed genes
confidently_expressed_genes_at1 <- list(
    coad = get_tissue_genes("coad", c("rectum", "colon"), min_expr = 1),
    luad = get_tissue_genes("luad", "lung", min_expr = 1),
    paad = get_tissue_genes("paad", "pancreas", min_expr = 1),
    skcm = get_tissue_genes("skcm", "skin", min_expr = 1),
    mm = get_tissue_genes("mm", c("blood", "bone marrow"), min_expr = 1)
)

confidently_expressed_genes_at0 <- list(
    coad = get_tissue_genes("coad", c("rectum", "colon"), min_expr = 0),
    luad = get_tissue_genes("luad", "lung", min_expr = 0),
    paad = get_tissue_genes("paad", "pancreas", min_expr = 0),
    skcm = get_tissue_genes("skcm", "skin", min_expr = 0),
    mm = get_tissue_genes("mm", c("blood", "bone marrow"), min_expr = 0)
)

genes_to_remove <- list(
    coad = # TODO:
)




# make list of unexpressed genes
colon_removed <- get_tissue_genes("coad", c("rectum", "colon"), min_expr = 0)
colon_removed <- colon_removed[!(colon_removed %in% colon_genes)]
length(colon_removed)
cat(c("not_expressed_genes", colon_removed),
    file = "removed_from_colon.txt",
    sep = "\n")
lung_removed <- get_tissue_genes("luad", "lung", min_expr = 0)
lung_removed <- lung_removed[!(lung_removed %in% lung_genes)]
length(lung_removed)
cat(c("not_expressed_genes", lung_removed),
    file = "removed_from_lung.txt",
    sep = "\n")
pancreas_removed <- get_tissue_genes("paad", "pancreas", min_expr = 0)
pancreas_removed <- pancreas_removed[!(pancreas_removed %in% pancreas_genes)]
length(pancreas_removed)
cat(c("not_expressed_genes", pancreas_removed),
    file = "removed_from_pancreas.txt",
    sep = "\n")
skin_removed <- get_tissue_genes("skcm", "skin", min_expr = 0)
skin_removed <- skin_removed[!(skin_removed %in% skin_genes)]
length(skin_removed)
cat(c("not_expressed_genes", skin_removed),
    file = "removed_from_skin.txt",
    sep = "\n")
blood_removed <- get_tissue_genes("mm", c("blood", "bone marrow"), min_expr = 0)
blood_removed <- blood_removed[!(blood_removed %in% blood_genes)]
length(blood_removed)
cat(c("not_expressed_genes", blood_removed),
    file = "removed_from_blood.txt",
    sep = "\n")

# it is likely better to removed genes from data that appear in the
# "removed_from_....txt" instead of only keeping genes in the "..._min1.txt"
# files so that genes are not removed because of imperfect gene naming


#### ---- Plotting ---- ####

# show number of genes per tissue
tibble(cancer = c("colon/rectum", "pancreas", "lung", "skin", "blood/bone marrow"),
       num_genes = c(length(colon_genes), length(pancreas_genes),
                     length(lung_genes), length(skin_genes),
                     length(blood_genes))) %>%
    ggplot(aes(x = cancer, y = num_genes)) +
    geom_col(aes(fill = cancer)) +
    theme_classic() +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.1))) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    labs (x = "", y = "number of genes expressed",
          title = "Gene expression per tissue")
ggsave("number_of_genes_barplot.jpeg",
       width = 4.5, height = 6.8, dpi = "retina")

# number of mutations removed from each cancer
mut_dat <- readRDS("../intermediate_files/processed_mutation_data.tib")

coad_removed <- mut_dat %>%
    filter(cancer == "coad") %>%
    filter(!(gene %in% colon_genes))
paad_removed <- mut_dat %>%
    filter(cancer == "paad") %>%
    filter(!(gene %in% pancreas_genes))
luad_removed <- mut_dat %>%
    filter(cancer == "luad") %>%
    filter(!(gene %in% lung_genes))
skcm_removed <- mut_dat %>%
    filter(cancer == "skcm") %>%
    filter(!(gene %in% skin_genes))
# TODO: add MM when it is in the mutation data

tibble(cancer = c("COAD", "PAAD", "LUAD", "SKCM"),
       num_muts = c(nrow(coad_removed),
                    nrow(paad_removed),
                    nrow(luad_removed),
                    nrow(skcm_removed)),
       num_genes = c(length(unique(coad_removed$gene)),
                     length(unique(paad_removed$gene)),
                     length(unique(luad_removed$gene)),
                     length(unique(skcm_removed$gene)))) %>%
    ggplot(aes(x = cancer, y = num_muts)) +
    geom_col(aes(fill = cancer)) +
    geom_text(aes(label = num_genes), nudge_y = 3000) +
    theme_classic() +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.1))) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    labs (x = "", y = "number of mutations removed",
          title = "Mutations removed per cancer",
          subtitle = "labeled with number of unique genes removed")
ggsave("removed_genes_barplot.jpeg",
       width = 4.5, height = 6.8, dpi = "retina")


coad_removed <- mut_dat %>%
    filter(cancer == "coad") %>%
    filter(!(gene %in% colon_genes)) %>%
    group_by(gene, cancer) %>%
    summarise(n_muts = n_distinct(sampleid)) %>%
    ungroup()
paad_removed <- mut_dat %>%
    filter(cancer == "paad") %>%
    filter(!(gene %in% pancreas_genes)) %>%
    group_by(gene, cancer) %>%
    summarise(n_muts = n_distinct(sampleid)) %>%
    ungroup()
luad_removed <- mut_dat %>%
    filter(cancer == "luad") %>%
    filter(!(gene %in% lung_genes)) %>%
    group_by(gene, cancer) %>%
    summarise(n_muts = n_distinct(sampleid)) %>%
    ungroup()
skcm_removed <- mut_dat %>%
    filter(cancer == "skcm") %>%
    filter(!(gene %in% skin_genes)) %>%
    group_by(gene, cancer) %>%
    summarise(n_muts = n_distinct(sampleid)) %>%
    ungroup()

removed_tib <- rbind(coad_removed, paad_removed) %>%
    rbind(luad_removed) %>%
    rbind(skcm_removed)

ggplot(removed_tib) +
    facet_grid(cancer ~ .) +
    geom_histogram(aes(n_muts, fill = cancer),
                   bins = 100,
                   position = "identity") +
    theme_classic() +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.1))) +
    scale_x_continuous(expand = expand_scale(mult = c(0, 0.1))) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    labs (x = "number of mutations", y = "count",
          title = "Number of mutations of removed genes")
ggsave("num_muts_removed_genes_hist.jpeg",
       width = 7, height = 7, dpi = "retina")



# histogram of mutational FREQUENCY of removed genes
coad_removed <- mut_dat %>%
    filter(cancer == "coad") %>%
    mutate(total_muts = n_distinct(sampleid)) %>%
    group_by(gene, cancer) %>%
    summarise(n_muts = n_distinct(sampleid) / unique(total_muts)) %>%
    ungroup() %>%
    filter(!(gene %in% colon_genes))
paad_removed <- mut_dat %>%
    filter(cancer == "paad") %>%
    mutate(total_muts = n_distinct(sampleid)) %>%
    group_by(gene, cancer) %>%
    summarise(n_muts = n_distinct(sampleid) / unique(total_muts)) %>%
    ungroup() %>%
    filter(!(gene %in% pancreas_genes))
luad_removed <- mut_dat %>%
    filter(cancer == "luad") %>%
    mutate(total_muts = n_distinct(sampleid)) %>%
    group_by(gene, cancer) %>%
    summarise(n_muts = n_distinct(sampleid) / unique(total_muts)) %>%
    ungroup() %>%
    filter(!(gene %in% lung_genes))
skcm_removed <- mut_dat %>%
    filter(cancer == "skcm") %>%
    mutate(total_muts = n_distinct(sampleid)) %>%
    group_by(gene, cancer) %>%
    summarise(n_muts = n_distinct(sampleid) / unique(total_muts)) %>%
    ungroup() %>%
    filter(!(gene %in% skin_genes))

removed_tib <- rbind(coad_removed, paad_removed) %>%
    rbind(luad_removed) %>%
    rbind(skcm_removed)

ggplot(removed_tib) +
    facet_grid(cancer ~ .) +
    geom_histogram(aes(n_muts, fill = cancer),
                   bins = 100,
                   position = "identity") +
    theme_classic() +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.1))) +
    scale_x_continuous(expand = expand_scale(mult = c(0, 0.1))) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +
    labs (x = "frequency of mutation", y = "count",
          title = "Frequency of mutations of removed genes")
ggsave("freq_muts_removed_genes_hist.jpeg",
       width = 7, height = 7, dpi = "retina")
