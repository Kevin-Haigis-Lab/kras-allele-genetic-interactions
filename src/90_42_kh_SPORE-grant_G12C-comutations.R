
TABLES_DIR <- "90_42_kh_SPORE-grant_G12C-comutations"
reset_table_directory(TABLES_DIR)

coad_comutation_interactions <- genetic_interaction_df %>%
    filter(cancer == "COAD")

g12c_incr_comutation <- coad_comutation_interactions %>%
    filter(allele == "G12C" & genetic_interaction == "comutation")

g12c_incr_comutation %>%
    select(hugo_symbol, allele, p_val, odds_ratio, n11, gene_freq) %>%
    arrange(-odds_ratio) %>%
    knitr::kable()
#> |hugo_symbol |allele |     p_val| odds_ratio| n11| gene_freq|
#> |:-----------|:------|---------:|----------:|---:|---------:|
#> |PHF8        |G12C   | 0.0052697|  11.869565|   3| 0.0101961|
#> |RBM27       |G12C   | 0.0052697|  11.869565|   3| 0.0101961|
#> |DCLK3       |G12C   | 0.0013752|  11.343798|   4| 0.0141176|
#> |HAO1        |G12C   | 0.0065735|  10.828607|   3| 0.0109804|
#> |MCM4        |G12C   | 0.0065735|  10.828607|   3| 0.0109804|
#> |PSG3        |G12C   | 0.0065735|  10.828607|   3| 0.0109804|
#> |PSG6        |G12C   | 0.0065735|  10.828607|   3| 0.0109804|
#> |MAGEE1      |G12C   | 0.0080536|   9.954203|   3| 0.0117647|
#> |SVIL        |G12C   | 0.0084843|   6.392449|   4| 0.0227451|
#> |ZNF462      |G12C   | 0.0034579|   6.282367|   5| 0.0290196|
#> |SCN1A       |G12C   | 0.0043776|   5.908361|   5| 0.0305882|
#> |ZNF804B     |G12C   | 0.0081810|   4.943838|   5| 0.0269773|
#> |MTOR        |G12C   | 0.0065102|   3.353245|   8| 0.0306574|

other_incr_comutation <- coad_comutation_interactions %>%
    filter(allele != "G12C" & genetic_interaction == "comutation")

g12c_incr_comutation %>%
    filter(!hugo_symbol %in% other_incr_comutation$hugo_symbol) %>%
    select(hugo_symbol, allele, p_val, odds_ratio, n11, gene_freq) %>%
    arrange(-odds_ratio) %>%
    knitr::kable()
#> |hugo_symbol |allele |     p_val| odds_ratio| n11| gene_freq|
#> |:-----------|:------|---------:|----------:|---:|---------:|
#> |PHF8        |G12C   | 0.0052697|  11.869565|   3| 0.0101961|
#> |RBM27       |G12C   | 0.0052697|  11.869565|   3| 0.0101961|
#> |DCLK3       |G12C   | 0.0013752|  11.343798|   4| 0.0141176|
#> |HAO1        |G12C   | 0.0065735|  10.828607|   3| 0.0109804|
#> |MCM4        |G12C   | 0.0065735|  10.828607|   3| 0.0109804|
#> |PSG3        |G12C   | 0.0065735|  10.828607|   3| 0.0109804|
#> |PSG6        |G12C   | 0.0065735|  10.828607|   3| 0.0109804|
#> |MAGEE1      |G12C   | 0.0080536|   9.954203|   3| 0.0117647|
#> |SVIL        |G12C   | 0.0084843|   6.392449|   4| 0.0227451|
#> |ZNF462      |G12C   | 0.0034579|   6.282367|   5| 0.0290196|
#> |SCN1A       |G12C   | 0.0043776|   5.908361|   5| 0.0305882|
#> |ZNF804B     |G12C   | 0.0081810|   4.943838|   5| 0.0269773|
#> |MTOR        |G12C   | 0.0065102|   3.353245|   8| 0.0306574|


#### ---- Top 5 most frequently comutated gene with KRAS G12C ---- ####

blacklist_genes <- c("KRAS", "TTN", "NEB")

num_g12c_mutants <- cancer_coding_muts_df %>%
    filter(ras_allele == "KRAS_G12C" & cancer == "COAD") %>%
    filter(target %in% c("exome", "genome")) %>%
    pull(tumor_sample_barcode) %>%
    n_distinct()

cancer_coding_muts_df %>%
    filter(ras_allele == "KRAS_G12C" & cancer == "COAD") %>%
    filter(target %in% c("exome", "genome")) %>%
    count(hugo_symbol, name = "num_comuts", sort = TRUE) %>%
    filter(!hugo_symbol %in% blacklist_genes) %>%
    filter(!str_detect(hugo_symbol, "^MUC|^HLA")) %>%
    top_n(10, wt = num_comuts) %>%
    mutate(percent_comut_with_g12c = round(num_comuts / !!num_g12c_mutants * 100)) %>%
    knitr::kable(format = "pandoc")




#### ---- G12C-specific comutants ---- ####

genetic_interaction_df %>%
    filter(cancer == "COAD" | cancer == "LUAD") %>%
    filter(genetic_interaction == "comutation") %>%
    group_by(cancer, hugo_symbol) %>%
    filter(n_distinct(kras_allele) == 1) %>%
    ungroup() %>%
    filter(kras_allele == "KRAS_G12C") %>%
    mutate(freq_in_g12c = n11 / (n11 + n01)) %>%
    select(cancer,
           gene = hugo_symbol,
           p_value = p_val,
           odds_ratio,
           gene_mutational_freq = gene_freq,
           freq_in_g12c_tumors = freq_in_g12c) %>%
    write_tsv(table_path(TABLES_DIR,
                         "g12c_sepcific_comutation_frequencies.tsv"))
