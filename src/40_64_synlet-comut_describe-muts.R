# Describe the mutations to TP53 in COAD lines and SMAD4 in PAAD lines.

srouce(file.path("src", "40_62_genes-for-synlet-comut.R"))


cell_lines_from_analysis <- synlet_comut_model_res %>%
    right_join(masking_hits, by = c("cancer", "allele", "hugo_symbol")) %>%
    group_by(cancer) %>%
    slice(1) %>%
    ungroup() %>%
    select(cancer, data) %>%
    unnest(data) %>%
    distinct(cancer, dep_map_id, kras_allele)

clines_mut_data <- ccle_mutations %>%
    inner_join(cell_lines_from_analysis, by = c("dep_map_id"))

genes_to_show_mutations <- masking_hits %>%
    distinct(cancer, other_gene) %>%
    rename(hugo_symbol = other_gene)

clines_mut_data %>%
    inner_join(genes_to_show_mutations, by = c("cancer", "hugo_symbol")) %>%
    select(cancer, hugo_symbol, variant_classification, protein_change,
           is_deleterious, is_tcga_hotspot, is_cosmic_hotspot) %>%
    knitr::kable(format = "markdown")

#> |cancer |hugo_symbol |variant_classification |protein_change |is_deleterious |is_tcga_hotspot |is_cosmic_hotspot |
#> |:------|:-----------|:----------------------|:--------------|:--------------|:---------------|:-----------------|
#> |PAAD   |SMAD4       |Missense_Mutation      |p.R100T        |FALSE          |FALSE           |TRUE              |
#> |PAAD   |SMAD4       |Nonsense_Mutation      |p.R135*        |TRUE           |TRUE            |TRUE              |
#> |PAAD   |SMAD4       |Splice_Site            |NA             |TRUE           |FALSE           |TRUE              |
#> |PAAD   |SMAD4       |Frame_Shift_Del        |p.QNGFT250fs   |TRUE           |FALSE           |TRUE              |
#> |PAAD   |SMAD4       |Frame_Shift_Del        |p.LD318fs      |TRUE           |FALSE           |TRUE              |
#> |PAAD   |SMAD4       |Missense_Mutation      |p.D319E        |FALSE          |FALSE           |TRUE              |
#> |PAAD   |SMAD4       |Frame_Shift_Del        |p.VHK330fs     |TRUE           |FALSE           |TRUE              |
#> |PAAD   |SMAD4       |Frame_Shift_Del        |p.P514fs       |TRUE           |FALSE           |TRUE              |
#> |COAD   |TP53        |Missense_Mutation      |p.P309S        |FALSE          |FALSE           |TRUE              |
#> |COAD   |TP53        |Missense_Mutation      |p.P278H        |FALSE          |TRUE            |TRUE              |
#> |COAD   |TP53        |Missense_Mutation      |p.R273L        |FALSE          |TRUE            |TRUE              |
#> |COAD   |TP53        |Missense_Mutation      |p.R273H        |FALSE          |TRUE            |TRUE              |
#> |COAD   |TP53        |Missense_Mutation      |p.R249S        |FALSE          |TRUE            |TRUE              |
#> |COAD   |TP53        |Missense_Mutation      |p.R248W        |FALSE          |TRUE            |TRUE              |
#> |COAD   |TP53        |Missense_Mutation      |NA             |FALSE          |TRUE            |FALSE             |
#> |COAD   |TP53        |Missense_Mutation      |p.R213L        |FALSE          |TRUE            |TRUE              |
#> |COAD   |TP53        |Nonsense_Mutation      |p.R213*        |TRUE           |TRUE            |TRUE              |
#> |COAD   |TP53        |Nonsense_Mutation      |p.E204*        |TRUE           |TRUE            |TRUE              |
#> |COAD   |TP53        |Missense_Mutation      |p.R196P        |FALSE          |TRUE            |TRUE              |
#> |COAD   |TP53        |Missense_Mutation      |p.P190L        |FALSE          |TRUE            |TRUE              |
#> |COAD   |TP53        |Splice_Site            |NA             |TRUE           |TRUE            |TRUE              |
#> |COAD   |TP53        |Missense_Mutation      |p.H179R        |FALSE          |TRUE            |TRUE              |
#> |COAD   |TP53        |Missense_Mutation      |p.R158L        |FALSE          |TRUE            |TRUE              |
#> |COAD   |TP53        |Frame_Shift_Del        |p.P72fs        |TRUE           |TRUE            |TRUE              |
#> |COAD   |TP53        |Nonsense_Mutation      |p.E51*         |TRUE           |TRUE            |TRUE              |
#> |COAD   |TP53        |Splice_Site            |NA             |TRUE           |FALSE           |FALSE             |
