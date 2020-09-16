# SETUP ----
sq_dna <- construct_sq_dna(
  c("ATCTTGAAG", "CATATGCGCTA", "ACGTGTCGA", ""), is_clean = TRUE)
sq_dna_compl <- construct_sq_dna(
  c("TAGAACTTC", "GTATACGCGAT", "TGCACAGCT", ""), is_clean = TRUE)
sq_rna <- construct_sq_rna(
  c("UAGUAACCGUAAGCG", "UAGUCC--UA-G"), is_clean = TRUE)
sq_rna_compl <- construct_sq_rna(
  c("AUCAUUGGCAUUCGC", "AUCAGG--AU-C"), is_clean = TRUE)

# PROTOTYPE PRESERVATION ----
test_that("complement() preserves all attributes of original vector", {
  expect_vector(complement(sq_dna),
                ptype = vec_ptype(sq_dna),
                size = vec_size(sq_dna))
  expect_vector(complement(sq_rna),
                ptype = vec_ptype(sq_rna),
                size = vec_size(sq_rna))
})

# VALUE COMPUTATION ----
# NOTE: used as.character() because hypothetically one value might have
#  multiple equivalent representations
test_that("complement() returns correct complement value for complement-only characters", {
  expect_equivalent(
    as.character(complement(sq_dna)),
    as.character(sq_dna_compl)
  )
})

test_that("complement() returns correct complement value for complement-less characters", {
  expect_equivalent(
    as.character(complement(sq_rna)),
    as.character(sq_rna_compl)
  )
})

# CANCELLING UPON DOUBLE USAGE ----
test_that("double use of complement() returns original value", {
  expect_identical(complement(complement(sq_dna)), sq_dna)
  expect_identical(complement(complement(sq_rna)), sq_rna)
})

# SHORTHAND FUNCTIONS ----
test_that("complement_dna() return identical value as complement() for DNA sequence", {
  expect_identical(complement(sq_dna), complement_dna(sq_dna))
})
test_that("complement_rna() return identical value as complement() for RNA sequence", {
  expect_identical(complement(sq_rna), complement_rna(sq_rna))
})
test_that("complement_dna() fail for RNA sequence", {
  expect_error(complement_dna(sq_rna))
})
test_that("complement_rna() fail for DNA sequence", {
  expect_error(complement_rna(sq_dna))
})
