
test_that("The KEGG parsing system works as before (snapshot).", {
  testthat::local_edition(3)
  p <- here::here("data", "kegg-pathways", "hsa04014_Ras-signaling-pathway.xml")
  g <- parse_and_annotate_kegg_kgml(p)
  expect_snapshot_value(as_tibble(g, active = "nodes"), style = "serialize")
  expect_snapshot_value(as_tibble(g, active = "edges"), style = "serialize")
})
