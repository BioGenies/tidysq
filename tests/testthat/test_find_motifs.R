sq_ami <- construct_sq(c("AGNTYIB", "MATEGILI", "MIPADHICA"), type = 'ami')
sq_nuc <- construct_sq(c("CTGAATGCAGT", "ATGCCGT", "CAGACCANN"), type = 'nuc')

test_that("find_motifs detects correctly motif that is single unambiguous amino acid in sequences", {
  expect_equal(dim(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "A"))[1],
               4)
  expect_equal(dim(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "C"))[1],
               1)
})

test_that("find_motifs detects correctly motif that is single unambiguous nucleotide in sequences", {
  expect_equal(dim(find_motifs(sq_nuc, c("sq1", "sq2", "sq3"), "A"))[1],
               7)
  expect_equal(dim(find_motifs(sq_nuc, c("sq1", "sq2", "sq3"), "C"))[1],
               7)
})

test_that("find_motifs detects correctly motif that is single ambiguous amino acid in sequences", {
  expect_equal(dim(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "B"))[1],
               3)
  expect_equal(dim(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "J"))[1],
               6)
})

test_that("find_motifs detects correctly motif that is single ambiguous nucleotide in sequences", {
  expect_equal(dim(find_motifs(sq_nuc, c("sq1", "sq2", "sq3"), "W"))[1],
               12)
  expect_equal(dim(find_motifs(sq_nuc, c("sq1", "sq2", "sq3"), "K"))[1],
               11)
})

test_that("find_motifs detects correctly leading letters of amino acid sequences using '^'", {
  expect_equal(dim(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "^A"))[1],
               1)
  expect_equal(dim(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "^M"))[1],
               2)
})

test_that("find_motifs detects correctly leading letters of nucleotide sequences using '^'", {
  expect_equal(dim(find_motifs(sq_nuc, c("sq1", "sq2", "sq3"), "^C"))[1],
               2)
  expect_equal(dim(find_motifs(sq_nuc, c("sq1", "sq2", "sq3"), "^A"))[1],
               1)
})

test_that("find_motifs detects correctly letters at the end of amino acid sequences using '$'", {
  expect_equal(dim(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "B$"))[1],
               1)
  expect_equal(dim(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "I$"))[1],
               1)
  expect_equal(dim(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "A$"))[1],
               1)
})

test_that("find_motifs detects correctly letters at the end of nucleotide sequences using '$'", {
  expect_equal(dim(find_motifs(sq_nuc, c("sq1", "sq2", "sq3"), "T$"))[1],
               2)
  expect_equal(dim(find_motifs(sq_nuc, c("sq1", "sq2", "sq3"), "N$"))[1],
               3)
})

test_that("find_motifs detects correctly multiple-letter motifs in amino acid sequences", {
  expect_equal(dim(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "TY"))[1],
               1)
  expect_equal(dim(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "PAD"))[1],
               1)
  expect_equal(dim(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "GILI"))[1],
               1)
  expect_equal(dim(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), "RARARARA"))[1],
               0)
})

test_that("find_motifs detects correctly multiple-letter motifs in nucleotide sequences", {
  expect_equal(dim(find_motifs(sq_nuc, c("sq1", "sq2", "sq3"), "CA"))[1],
               3)
  expect_equal(dim(find_motifs(sq_nuc, c("sq1", "sq2", "sq3"), "CAN"))[1],
               3)
  expect_equal(dim(find_motifs(sq_nuc, c("sq1", "sq2", "sq3"), "GAAT"))[1],
               1)
  expect_equal(dim(find_motifs(sq_nuc, c("sq1", "sq2", "sq3"), "RARARARA"))[1],
               0)
})

test_that("find_motifs detects correctly multiple motifs in amino acid sequences", {
  expect_equal(dim(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), c("TY", "TYM")))[1],
               1)
  expect_equal(dim(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), c("PAD", "TY", "GILI")))[1],
               3)
  expect_equal(dim(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), c("GILI", "PAD")))[1],
               2)
  expect_equal(dim(find_motifs(sq_ami, c("sq1", "sq2", "sq3"), c("RARARARA", "GILI")))[1],
               1)
})

test_that("find_motifs detects correctly multiple motifs in nucleotide sequences", {
  expect_equal(dim(find_motifs(sq_nuc, c("sq1", "sq2", "sq3"), c("CA", "CAN")))[1],
               6)
  expect_equal(dim(find_motifs(sq_nuc, c("sq1", "sq2", "sq3"), c("CAN", "BRRMAA")))[1],
               3)
  expect_equal(dim(find_motifs(sq_nuc, c("sq1", "sq2", "sq3"), c("GAAT", "CAN")))[1],
               4)
  expect_equal(dim(find_motifs(sq_nuc, c("sq1", "sq2", "sq3"), c("RARARARA", "BYBYBYBYBYBY")))[1],
               0)
})