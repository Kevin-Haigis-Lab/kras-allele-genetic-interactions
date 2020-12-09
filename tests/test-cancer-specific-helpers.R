
cancer_names <- c("COAD", "LUAD", "PAAD", "MM", "SKCM")

# test `get_hugo_from_depmap_ids()`
test_that("Cancer specific named vectors contain valid cancers.", {

  compare_against_cancer_names <- function(x) {
    expect_true(all(names(x) %in% cancer_names))
  }

  compare_against_cancer_names(cancer_oncogenes)
  compare_against_cancer_names(FRACTION_OF_TISSUE_THAT_ARE_CANCER)
  compare_against_cancer_names(cgc_cancer_regex)
})


test_that("The COSMIC CGC genes can be retrieved.", {
  expect_error(get_cgc_genes("NOTCANCER"))
  for (cancer in cancer_names[cancer_names != "SKCM"]) {
    genes <- get_cgc_genes(cancer)
    expect_true(all(is.character(genes)))
    expect_true(!any(is.na(genes)))
  }
})