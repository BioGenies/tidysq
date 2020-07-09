sq_ami <- construct_sq(c("AGNTYIB", "MATEGILI", "MIPADHICA"), type = 'ami')
sq_dna <- construct_sq(c("CTGAATGCAGT", "ATGCCGT", "CAGACCANN"), type = 'dna')
sq_rna <- construct_sq(c("UUACGACUU", "UUAAGCGC", "ACUAAGACCA"), type = 'rna')

test_that("find_motifs detects correctly motif that is single unambiguous amino acid in sequences", {
  expect_equal(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "A")[["start"]],
               c(1, 2, 4, 9))
  expect_equal(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "C")[["start"]],
               c(start = 8))
})

test_that("find_motifs detects correctly motif that is single unambiguous DNA nucleotide in sequences", {
  expect_equal(find_motifs(sq_dna, c("sq1", "sq2", "sq3"), "A")[["start"]],
               c(4, 5, 9, 1, 2, 4, 7))
  expect_equal(find_motifs(sq_dna, c("sq1", "sq2", "sq3"), "C")[["start"]],
               c(1, 8, 4, 5, 1, 5, 6))
})

test_that("find_motifs detects correctly motif that is single unambiguous RNA nucleotide in sequences", {
  expect_equal(find_motifs(sq_rna, c("sq1", "sq2", "sq3"), "A")[["start"]],
               c(3, 6, 3, 4, 1, 4, 5, 7, 10))
  expect_equal(find_motifs(sq_rna, c("sq1", "sq2", "sq3"), "C")[["start"]],
               c(4, 7, 6, 8, 2, 8, 9))
})

test_that("find_motifs detects correctly motif that is single ambiguous amino acid in sequences", {
  expect_equal(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "B")[["start"]],
               c(3, 7, 5))
  expect_equal(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "J")[["start"]],
               c(6, 6, 7, 8, 2, 7))
})

test_that("find_motifs detects correctly motif that is single ambiguous DNA nucleotide in sequences", {
  expect_equal(find_motifs(sq_dna, c("sq1", "sq2", "sq3"), "W")[["start"]],
               c(2, 4, 5, 6, 9, 11, 1, 2, 7, 2, 4, 7))
  expect_equal(find_motifs(sq_dna, c("sq1", "sq2", "sq3"), "K")[["start"]],
               c(2, 3, 6, 7, 10, 11, 2, 3, 6, 7, 3))
})

test_that("find_motifs detects correctly motif that is single ambiguous RNA nucleotide in sequences", {
  expect_equal(find_motifs(sq_rna, c("sq1", "sq2", "sq3"), "B")[["start"]],
               c(1, 2, 4, 5, 7, 8, 9, 1, 2, 5, 6, 7, 8, 2, 3, 6, 8, 9))
  expect_equal(find_motifs(sq_rna, c("sq1", "sq2", "sq3"), "Y")[["start"]],
               c(1, 2, 4, 7, 8, 9, 1, 2, 6, 8, 2, 3, 8, 9))
})

test_that("find_motifs detects correctly leading letters of amino acid sequences using '^'", {
  expect_equal(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "^A")[["start"]],
               c(start = 1))
  expect_equal(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "^M")[["start"]],
               c(1, 1))
})

test_that("find_motifs detects correctly leading letters of DNA nucleotide sequences using '^'", {
  expect_equal(find_motifs(sq_dna, c("sq1", "sq2", "sq3"), "^C")[["start"]],
               c(1, 1))
  expect_equal(find_motifs(sq_dna, c("sq1", "sq2", "sq3"), "^A")[["start"]],
               c(start = 1))
})

test_that("find_motifs detects correctly leading letters of RNA nucleotide sequences using '^'", {
  expect_equal(find_motifs(sq_rna, c("sq1", "sq2", "sq3"), "^C")[["start"]],
               numeric(0))
  expect_equal(find_motifs(sq_rna, c("sq1", "sq2", "sq3"), "^A")[["start"]],
               c(start = 1))
  expect_equal(find_motifs(sq_rna, c("sq1", "sq2", "sq3"), "^UU")[["start"]],
               c(1, 1))
})

test_that("find_motifs detects correctly letters at the end of amino acid sequences using '$'", {
  expect_equal(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "B$")[["start"]],
               c(start = 7))
  expect_equal(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "I$")[["start"]],
               c(start = 8))
  expect_equal(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "A$")[["start"]],
               c(start = 9))
})

test_that("find_motifs detects correctly letters at the end of DNA nucleotide sequences using '$'", {
  expect_equal(find_motifs(sq_dna, c("sq1", "sq2", "sq3"), "T$")[["start"]],
               c(11, 7))
  expect_equal(find_motifs(sq_dna, c("sq1", "sq2", "sq3"), "N$")[["start"]],
               c(11, 7, 9))
})

test_that("find_motifs detects correctly letters at the end of RNA nucleotide sequences using '$'", {
  expect_equal(find_motifs(sq_rna, c("sq1", "sq2", "sq3"), "D$")[["start"]],
               c(9, 10))
  expect_equal(find_motifs(sq_rna, c("sq1", "sq2", "sq3"), "C$")[["start"]],
               c(start = 8))
})

test_that("find_motifs detects correctly multiple-letter motifs in amino acid sequences", {
  expect_equal(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "TY")[["start"]],
               c(start = 4))
  expect_equal(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "PAD")[["start"]],
               c(start = 3))
  expect_equal(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "GILI")[["start"]],
               c(start = 5))
  expect_equal(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "RARARARA")[["start"]],
               numeric(0))
})

test_that("find_motifs detects correctly multiple-letter motifs in nucleotide sequences", {
  expect_equal(find_motifs(sq_dna, c("sq1", "sq2", "sq3"), "CA")[["start"]],
               c(8, 1, 6))
  expect_equal(find_motifs(sq_rna, c("sq1", "sq2", "sq3"), "NCA")[["start"]],
               c(start = 8))
  expect_equal(find_motifs(sq_dna, c("sq1", "sq2", "sq3"), "GAAT")[["start"]],
               c(start = 3))
  expect_equal(find_motifs(sq_rna, c("sq1", "sq2", "sq3"), "RARARARA")[["start"]],
               numeric(0))
})

test_that("find_motifs detects correctly multiple motifs in amino acid sequences", {
  expect_equal(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), c("TY", "TYM"))[["start"]],
               c(start = 4))
  expect_equal(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), c("PAD", "TY", "GILI"))[["start"]],
               c(3, 4, 5))
  expect_equal(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), c("GILI", "PAD"))[["start"]],
               c(5, 3))
  expect_equal(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), c("RARARARA", "GILI"))[["start"]],
               c(start = 5))
})

test_that("find_motifs detects correctly multiple motifs in nucleotide sequences", {
  expect_equal(find_motifs(sq_dna, c("sq1", "sq2", "sq3"), c("CA", "CAN"))[["start"]],
               c(8, 1, 6, 8, 1, 6))
  expect_equal(find_motifs(sq_rna, c("sq1", "sq2", "sq3"), c("CAN", "BRA"))[["start"]],
               c(4, 2, 3))
  expect_equal(find_motifs(sq_dna, c("sq1", "sq2", "sq3"), c("GAAT", "CAN"))[["start"]],
               c(3, 8, 1, 6))
  expect_equal(find_motifs(sq_rna, c("sq1", "sq2", "sq3"), c("RARARARA", "BYBYBYBYBYBY"))[["start"]],
               numeric(0))
})