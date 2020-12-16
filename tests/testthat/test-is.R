# SETUP ----
sq_dna_bsc <- sq(c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT"),
                 alphabet = "dna_bsc")
sq_dna_ext <- sq(c("NARTYVTCY", "", "ATKCYGDD", "", "DNAKYTD"),
                 alphabet = "dna_ext")
sq_rna_bsc <- sq(c("UAUCAGU-A-GU-CA", "CUG-A-CUGAG-CC", "-CUG-AGAGUA-"),
                 alphabet = "rna_bsc")
sq_rna_ext <- sq(c("", "KBS-UVW-AWWWG", "YGHHH-", "-CRASH", "MND-KUUBV-MY-"),
                 alphabet = "rna_ext")
sq_ami_bsc <- sq(c("ACEH", "PASAI", "MALACCA", "SIAK"),
                 alphabet = "ami_bsc")
sq_ami_ext <- sq(c("OUTLANDISH", "UNSTRUCTURIZED", "FEAR"),
                 alphabet = "ami_ext")
sq_unt <- sq(c("VIP01", "VIP02", "VIP04", "MISSING_ONE"),
             alphabet = "unt")
sq_atp <- sq(c("mAmYmY", "nbAnsAmA", ""),
             alphabet = c("mA", "mY", "nbA", "nsA"))

# SQ CHECKING ----
test_that("any sq type is considered to be of sq class", {
  expect_true(is.sq(sq_dna_bsc))
  expect_true(is.sq(sq_dna_ext))
  expect_true(is.sq(sq_rna_bsc))
  expect_true(is.sq(sq_rna_ext))
  expect_true(is.sq(sq_ami_bsc))
  expect_true(is.sq(sq_ami_ext))
  expect_true(is.sq(sq_unt))
  expect_true(is.sq(sq_atp))
})

test_that("other objects are not considered sq", {
  expect_false(is.sq(c(1,2,3)))
  expect_false(is.sq(LETTERS))
  expect_false(is.sq(TRUE))
  expect_false(is.sq(NULL))
  # That's a function, not an sq object!
  expect_false(is.sq(sq))
})

# EXACT CLASS CHECKING ----
test_that("objects can be checked for exact type", {
  expect_true(is.sq_dna_bsc(sq_dna_bsc))
  expect_false(is.sq_dna_bsc(sq_dna_ext))
  expect_true(is.sq_dna_ext(sq_dna_ext))
  expect_true(is.sq_rna_bsc(sq_rna_bsc))
  expect_true(is.sq_rna_ext(sq_rna_ext))
  expect_false(is.sq_rna_ext(sq_unt))
  expect_true(is.sq_ami_bsc(sq_ami_bsc))
  expect_true(is.sq_ami_ext(sq_ami_ext))
  expect_true(is.sq_unt(sq_unt))
  expect_false(is.sq_unt(sq_atp))
  expect_true(is.sq_atp(sq_atp))
  expect_false(is.sq_atp(sq_unt))
})

# GENERALIZED CLASS CHECKING ----
test_that("DNA, RNA and amino acid sequences can be checked regardless of basic or extended type", {
  expect_true(is.sq_dna(sq_dna_bsc))
  expect_true(is.sq_dna(sq_dna_ext))
  expect_true(is.sq_rna(sq_rna_bsc))
  expect_false(is.sq_rna(sq_dna_bsc))
  expect_false(is.sq_rna(sq_unt))
  expect_true(is.sq_ami(sq_ami_ext))
  expect_false(is.sq_ami(sq_dna_bsc))
})
