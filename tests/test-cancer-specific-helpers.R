
cancer_names <- c("COAD", "LUAD", "PAAD", "MM", "SKCM")

# test `get_hugo_from_depmap_ids()`
test_that("Cancer specific named vectors contain valid cancers.", {
  expect_true(all(names(cancer_oncogenes) %in% cancer_names))
  expect_true(all(names(FRACTION_OF_TISSUE_THAT_ARE_CANCER) %in% cancer_names))
  expect_true(all(names(cgc_cancer_regex) %in% cancer_names))
})


test_that("The COSMIC CGC genes can be retrieved.", {
  expect_error(get_cgc_genes("NOTCANCER"))
  for (cancer in cancer_names[cancer_names != "SKCM"]) {
    genes <- get_cgc_genes(cancer)
    expect_true(all(is.character(genes)))
    expect_true(!any(is.na(genes)))
    # expect_snapshot_value(get_cgc_genes(cancer))
  }
})