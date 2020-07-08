
GRAPHS_DIR <- "90_42_kh_SPORE-grant_G12C-comutations"
reset_graph_directory(GRAPHS_DIR)

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