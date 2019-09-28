
# test `get_hugo_from_depmap_ids()`
test_that("Hugo symbols are extracted.", {
    test_genes <- c(
        "ABCG8 (64241)",
        "ACCSL (390110)",
        "ACTL7A (10881)",
        "ACTL7B (10880)",
        "ACTL9 (284382)",
        "ADAD1 (132612)",
        "ADAM18 (8749)",
        "ADAM2 (2515)",
        "ADAM20 (8748)",
        "ADAM30 (11085)"
    )

    expected_output <- c(
        "ABCG8",
        "ACCSL",
        "ACTL7A",
        "ACTL7B",
        "ACTL9",
        "ADAD1",
        "ADAM18",
        "ADAM2",
        "ADAM20",
        "ADAM30"
    )

    expect_length(get_hugo_from_depmap_ids(test_genes), length(test_genes))
    expect_equal(get_hugo_from_depmap_ids(test_genes), expected_output)
    expect_length(get_hugo_from_depmap_ids(NULL), 0)
})

# test `get_entrez_from_depmap_ids()`
test_that("Entrez IDs are extracted.", {
    test_genes <- c(
        "ABCG8 (64241)",
        "ACCSL (390110)",
        "ACTL7A (10881)",
        "ACTL7B (10880)",
        "ACTL9 (284382)",
        "ADAD1 (132612)",
        "ADAM18 (8749)",
        "ADAM2 (2515)",
        "ADAM20 (8748)",
        "ADAM30 (11085)"
    )

    expected_output <- c(
        "64241",
        "390110",
        "10881",
        "10880",
        "284382",
        "132612",
        "8749",
        "2515",
        "8748",
        "11085"
    )

    expect_length(get_entrez_from_depmap_ids(test_genes), length(test_genes))
    expect_equal(get_entrez_from_depmap_ids(test_genes), expected_output)
    expect_equal(get_entrez_from_depmap_ids(test_genes, convert_to_num = TRUE),
                 as.numeric(expected_output))
    expect_length(get_entrez_from_depmap_ids(NULL), 0)
})
