
testthat::test_that("Terms are properly italicized", {
  expect_equal(italicize_genes_in_terms("kras_allele"), "kras_allele")
  expect_equal(italicize_genes_in_terms("GENE"), "*GENE*")
  expect_equal(italicize_genes_in_terms(c("GENE", "kras_allele")), c("*GENE*", "kras_allele"))
  expect_equal(italicize_genes_in_terms("GENE1:GENE2"), "*GENE1*:*GENE2*")
})
