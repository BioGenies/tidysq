# SETUP ----
sq_ami <- sq(c("", "MATEGILI", "MIPADHICA"), alphabet = "ami_ext")
sq_dna <- sq(c("", "ATGCCGT", "T", "", "TTCGATCAGGC"), alphabet = "dna_bsc")
sq_rna <- sq(c("UUCAGC", "UAGUACCGA", "CAGGGGGGA"), alphabet = "rna_bsc")
sq_empty <- sq(character(), alphabet = "dna_ext")

# RETURNING LOGICAL VECTOR ----
test_that("is_empty_sq() returns logical vector", {
  expect_vector(is_empty_sq(sq_ami),
                ptype = logical(),
                size = length(sq_ami))
  expect_vector(is_empty_sq(sq_dna),
                ptype = logical(),
                size = length(sq_dna))
  expect_vector(is_empty_sq(sq_rna),
                ptype = logical(),
                size = length(sq_rna))
  expect_vector(is_empty_sq(sq_empty),
                ptype = logical(),
                size = length(sq_empty))
})

# ERROR FOR NON-SQ OBJECTS ----
test_sq_only(is_empty_sq)

# VALUE COMPUTATION
test_that("is_empty_sq() correctly computes return value", {
  expect_equal(is_empty_sq(sq_ami), c(TRUE, FALSE, FALSE))
  expect_equal(is_empty_sq(sq_dna), c(TRUE, FALSE, FALSE, TRUE, FALSE))
  expect_equal(is_empty_sq(sq_rna), c(FALSE, FALSE, FALSE))
  expect_equal(is_empty_sq(sq_empty), logical())
})
