# SETUP ----
sq_ami <- construct_sq(c("", "MATEGILI", "MIPADHICA"), type = 'ami')
sq_dna <- construct_sq(c("", "ATGCCGT", "T", "", "TTCGATCAGGC"), type = 'dna')
sq_rna <- construct_sq(c("UUCAGC", "UAGUACCGA", "CAGGGGGGA"), type = "rna")
sq_empty <- construct_sq(character(), type = "dna")

# RETURNING LOGICAL VECTOR ----
test_that("is_null_sq() returns logical vector", {
  expect_vector(is_null_sq(sq_ami),
                ptype = logical(),
                size = length(sq_ami))
  expect_vector(is_null_sq(sq_dna),
                ptype = logical(),
                size = length(sq_dna))
  expect_vector(is_null_sq(sq_rna),
                ptype = logical(),
                size = length(sq_rna))
  expect_vector(is_null_sq(sq_empty),
                ptype = logical(),
                size = length(sq_empty))
})

# VALUE COMPUTATION
test_that("is_null_sq() correctly computes return value", {
  expect_equivalent(is_null_sq(sq_ami), c(TRUE, FALSE, FALSE))
  expect_equivalent(is_null_sq(sq_dna), c(TRUE, FALSE, FALSE, TRUE, FALSE))
  expect_equivalent(is_null_sq(sq_rna), c(FALSE, FALSE, FALSE))
  expect_equivalent(is_null_sq(sq_empty), logical())
})
