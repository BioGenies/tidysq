sq_ami <- construct_sq(c("AGNTYIKFGGAYTIB", "MATEGILIAADGYTWIL", "MIPADHICAANGIENAGIK"), type = 'ami')
sq_dna <- construct_sq(c("CTGAATGCAGTACCGTAAT", "ATGCCGTAAATGCCAT", "CAGACCANNNATAG"), type = 'dna')
sq_rna <- construct_sq(c("GGCUGCGGGACUGAGGC", "UUCAUGCGGCUAGGGCU", "UAGGCGAGCAGAGUAG"), type = 'rna')

test_that("%has% detects correctly motif that is single unambiguous amino acid in sequences", {
  expect_equal(sq_ami %has% "A",
               c(TRUE, TRUE, TRUE))
  expect_equal(sq_ami %has% "C",
               c(FALSE, FALSE, TRUE))
})

test_that("%has% detects correctly motif that is single unambiguous DNA nucleotide in sequences", {
  expect_equal(sq_dna %has% "A",
               c(TRUE, TRUE, TRUE))
  expect_equal(sq_dna %has% "C",
               c(TRUE, TRUE, TRUE))
})

test_that("%has% detects correctly motif that is single unambiguous RNA nucleotide in sequences", {
  expect_equal(sq_rna %has% "A",
               c(TRUE, TRUE, TRUE))
  expect_equal(sq_rna %has% "C",
               c(TRUE, TRUE, TRUE))
})

test_that("%has% detects correctly motif that is single ambiguous amino acid in sequences", {
  expect_equal(sq_ami %has% "B", 
               c(TRUE, TRUE, TRUE))
  expect_equal(sq_ami %has% "J", 
               c(TRUE, TRUE, TRUE))
})

test_that("%has% detects correctly motif that is single ambiguous DNA nucleotide in sequences", {
  expect_equal(sq_dna %has% "W",
               c(TRUE, TRUE, TRUE))
  expect_equal(sq_dna %has% "K",
               c(TRUE, TRUE, TRUE))
})

test_that("%has% detects correctly motif that is single ambiguous RNA nucleotide in sequences", {
  expect_equal(sq_rna %has% "W",
               c(TRUE, TRUE, TRUE))
  expect_equal(sq_rna %has% "K",
               c(TRUE, TRUE, TRUE))
})

test_that("%has% detects correctly leading letters of amino acid sequences using '^'", {
  expect_equal(sq_ami %has% "^A",
               c(TRUE, FALSE, FALSE))
  expect_equal(sq_ami %has% "^M",
               c(FALSE, TRUE, TRUE))
})

test_that("%has% detects correctly leading letters of DNA nucleotide sequences using '^'", {
  expect_equal(sq_dna %has% "^C", 
               c(TRUE, FALSE, TRUE))
  expect_equal(sq_dna %has% "^A", 
               c(FALSE, TRUE, FALSE))
})

test_that("%has% detects correctly leading letters of RNA nucleotide sequences using '^'", {
  expect_equal(sq_rna %has% "^C", 
               c(FALSE, FALSE, FALSE))
  expect_equal(sq_rna %has% "^U", 
               c(FALSE, TRUE, TRUE))
})

test_that("%has% detects correctly letters at the end of amino acid sequences using '$'", {
  expect_equal(sq_ami %has% "K$", 
               c(FALSE, FALSE, TRUE))
  expect_equal(sq_ami %has% "L$", 
               c(FALSE, TRUE, FALSE))
})

test_that("%has% detects correctly letters at the end of DNA nucleotide sequences using '$'", {
  expect_equal(sq_dna %has% "T$",
               c(TRUE, TRUE, FALSE))
  expect_equal(sq_dna %has% "G$",
               c(FALSE, FALSE, TRUE))
})

test_that("%has% detects correctly letters at the end of RNA nucleotide sequences using '$'", {
  expect_equal(sq_rna %has% "Y$",
               c(TRUE, TRUE, FALSE))
  expect_equal(sq_rna %has% "G$",
               c(FALSE, FALSE, TRUE))
})

test_that("%has% detects correctly multiple-letter motifs in amino acid sequences", {
  expect_equal(sq_ami %has% "TY",
               c(TRUE, FALSE, FALSE))
  expect_equal(sq_ami %has% "GAY",
               c(TRUE, FALSE, FALSE))
  expect_equal(sq_ami %has% "MATE",
               c(FALSE, TRUE, FALSE))
})

test_that("%has% detects correctly multiple-letter motifs in DNA nucleotide sequences", {
  expect_equal(sq_dna %has% "CC",
               c(TRUE, TRUE, TRUE))
  expect_equal(sq_dna %has% "TAA",
               c(TRUE, TRUE, FALSE))
  expect_equal(sq_dna %has% "TAAT",
               c(TRUE, FALSE, FALSE))
})

test_that("%has% detects correctly multiple-letter motifs in RNA nucleotide sequences", {
  expect_equal(sq_rna %has% "CC",
               c(FALSE, FALSE, FALSE))
  expect_equal(sq_rna %has% "UAG",
               c(FALSE, TRUE, TRUE))
  expect_equal(sq_rna %has% "GGGA",
               c(TRUE, FALSE, FALSE))
})

test_that("%has% detects correctly multiple motifs in amino acid sequences", {
  expect_equal(sq_ami %has% c("TY", "GA"),
               c(TRUE, FALSE, FALSE))
  expect_equal(sq_ami %has% c("MATE", "TY"),
               c(FALSE, FALSE, FALSE))
  expect_equal(sq_ami %has% c("GI", "MATE"),
               c(FALSE, TRUE, FALSE))
})

test_that("%has% detects correctly multiple motifs in DNA nucleotide sequences", {
  expect_equal(sq_dna %has% c("CC", "AT"),
               c(TRUE, TRUE, TRUE))
  expect_equal(sq_dna %has% c("CC", "TAAT"),
               c(TRUE, FALSE, FALSE))
  expect_equal(sq_dna %has% c("CC", "TAA"),
               c(TRUE, TRUE, FALSE))
})

test_that("%has% detects correctly multiple motifs in RNA nucleotide sequences", {
  expect_equal(sq_rna %has% c("GC", "GA"),
               c(TRUE, FALSE, TRUE))
  expect_equal(sq_rna %has% c("UA", "CA"),
               c(FALSE, TRUE, TRUE))
  expect_equal(sq_rna %has% c("AUC", "ACC"),
               c(FALSE, FALSE, FALSE))
})

test_that("%has% works on a single sequence", {
  single_seq <- construct_sq("AGNTYIKFGGAYTIB", type = 'ami')
  expect_true(single_seq %has% "AG")
  expect_false(single_seq %has% "FBI")
})

