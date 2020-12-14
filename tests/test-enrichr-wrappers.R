
test_that("The EnrichR data bases are collected from the API.", {
  expect_true(nrow(all_enrichr_dbs) > 1)
  expect_true(ncol(all_enrichr_dbs) == 5)
})

test_that("The overlap column is split correctly.", {
  expect_equal(get_enrichr_overlap_int("12/34"), 12)
  expect_equal(get_enrichr_overlap_int(" 12/34"), 12)
})

test_that("Character list of genes is split correctly.", {
  expect_equal(enrichr_genes("ABC;DEF;GHI"), c("ABC", "DEF", "GHI"))
  expect_equal(enrichr_genes("ABC; DEF; GHI"), c("ABC", "DEF", "GHI"))
  expect_equal(enrichr_genes("ABC;"), c("ABC"))
  expect_equal(enrichr_genes("ABC"), c("ABC"))
  expect_equal(enrichr_genes("A"), c("A"))
  expect_equal(length(enrichr_genes("")), 0)
})
