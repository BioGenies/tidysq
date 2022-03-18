# SETUP ----
str_dna <- c("TTCAGGGCTAG", "CGATTGC", "CAGTTTA")
str_unt <- c("!!NV!!XFD!", "P!OQ!!-FI", "SP!!I-F-XXS!")
str_atp <- c("mAmYmY", "nbAnsAmA", "")

alph_atp <- c("mA", "mY", "nbA", "nsA")

sq_dna <- sq(str_dna, alphabet = "dna_bsc")
sq_unt <- sq(str_unt, alphabet = "unt", NA_letter = "!")
sq_atp <- sq(str_atp, alphabet = alph_atp)
sq_empty <- sq(character(), alphabet = "rna_bsc")
sq_one <- sq("OFUFQBODBUDAFNDBZGFG", alphabet = "ami_ext")

# PROTOTYPE PRESERVATION ----
test_that("collapse() preserves type and alphabet of original vector", {
  expect_vector(collapse(sq_dna),
                ptype = vec_ptype(sq_dna))
  expect_vector(collapse(sq_unt),
                ptype = vec_ptype(sq_unt))
  expect_vector(collapse(sq_atp),
                ptype = vec_ptype(sq_atp))
  expect_vector(collapse(sq_empty),
                ptype = vec_ptype(sq_empty))
})

test_that("collapse() returns sequence vector of length 1", {
  expect_vector(collapse(sq_dna),
                size = 1)
  expect_vector(collapse(sq_unt),
                size = 1)
  expect_vector(collapse(sq_atp),
                size = 1)
  expect_vector(collapse(sq_empty),
                size = 1)
})

# ERROR FOR NON-SQ OBJECTS ----
test_sq_only(collapse)

# VALUE COMPUTATION ----
test_that("collapse() returns sequence equal to all collapsed sequences", {
  withr::local_options(list(tidysq_NA_letter = "!"))
  expect_equal(as.character(collapse(sq_dna)),
               paste0(str_dna, collapse = ""))
  expect_equal(as.character(collapse(sq_unt)),
               paste0(str_unt, collapse = ""))
  expect_equal(as.character(collapse(sq_atp)),
               paste0(str_atp, collapse = ""))
  # Note that returned value isn't an empty vector, but an empty string
  expect_equal(as.character(collapse(sq_empty)),
               "")
})

# NO CHANGES FOR SQ OF LENGTH 1 ----
test_that("collapse() does nothing to sq objects of length 1", {
  expect_identical(collapse(sq_one),
                   sq_one)
})
