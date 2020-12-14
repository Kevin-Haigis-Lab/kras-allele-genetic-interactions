
test_that("The DNA nucleotides are correct.", {
  correct_nucs <- c("A", "C", "G", "T")
  expect_true(all(DNA_NUCLEOTIDES %in% correct_nucs))
  expect_true(all(correct_nucs %in% DNA_NUCLEOTIDES))
  expect_length(setdiff(DNA_NUCLEOTIDES, correct_nucs), 0)
  expect_length(setdiff(correct_nucs, DNA_NUCLEOTIDES), 0)
  expect_length(DNA_NUCLEOTIDES, 4)
  expect_equal(str_to_upper(DNA_NUCLEOTIDES), DNA_NUCLEOTIDES)
})

test_that("The correct nucleotide is replaced.", {
  dna <- Biostrings::DNAString("ACGTACGT")
  expect_equal(
    replace_nucleotide(dna, 2, "A"),
    Biostrings::DNAString("AAGTACGT")
  )
})

test_that("Codon translation.", {
  expect_equal(translate_codon(Biostrings::DNAString("ACT")), "T")
  expect_equal(translate_codon(Biostrings::DNAString("CCC")), "P")
  expect_equal(translate_codon(Biostrings::DNAString("AAA")), "K")
  expect_equal(translate_codon(Biostrings::DNAString("CGA")), "R")
  expect_equal(translate_codon(Biostrings::DNAString("CGG")), "R")
})

test_that("The mutation data frame is created properly.", {
  expect_output(enumerate_sbs_mutations("ATC"))
  expect_error(enumerate_sbs_mutations("ATCA"), "length")
  df <- enumerate_sbs_mutations("ATC")
  expect_equal(nrow(df), 13)
  expect_equal(ncol(df), 4)
  expect_equal(colnames(df), c("nuc", "pos", "codon", "amino_acid"))
  expect_true(all(df$amino_acid %in% LETTERS))
})
