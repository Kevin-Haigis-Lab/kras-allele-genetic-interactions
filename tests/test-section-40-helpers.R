
test_that("Data frames are returned by getters", {
  expect_s3_class(get_synthetic_lethal_data("COAD", "G12D"), "tbl_df")
  expect_s3_class(get_genetic_interaction_data("COAD", "G12D"), "tbl_df")
  expect_s3_class(get_overlapped_df("COAD", "G12D"), "tbl_df")

  gr <- get_overlapped_gr("COAD", "G12D")
  expect_s3_class(gr$graph, "tbl_graph")
  expect_s3_class(gr$data, "tbl_df")
})