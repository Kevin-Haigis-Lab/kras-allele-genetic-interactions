

cache("tcga_purity_ploidy", {
  data_file <- file.path(
    "data", "tcga",
    "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"
  )
  tcga_purity_ploidy <- read_tsv(data_file) %>%
    janitor::clean_names() %>%
    mutate(tumor_sample_barcode = str_remove(array, "-[:digit:]{2}$")) %>%
    select(-coverage_for_80_percent_power, -solution, -sample) %>%
    rename(full_id = array)
  rename(tcga_purity_ploidy)
})
