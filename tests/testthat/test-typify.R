# SETUP ----
sq_unt <- sq(c("PQNVIXFD", "PDOQXN-FI", "SPBI-F-XXS"), alphabet = "unt")
sq_atp <- sq(c("mAmYmY", "nbAnsAmA", ""),
             alphabet = c("mA", "mY", "nbA", "nsA"))
sq_dna <- sq(c("TTCAGGGCTAG", "CGATTGC", "CAGTTTA"), alphabet = "dna_bsc")
sq_empty <- sq(character(), alphabet = "dna_bsc")

# CORRECT PROTOTYPE OF RETURNED VALUE ----
test_that("typify() returns an sq object of desired type", {
  expect_vector(typify(sq_unt, "ami_ext"),
                ptype = sq_ptype(CPP_get_standard_alphabet("ami_ext"), "ami_ext"),
                size = vec_size(sq_unt))
  expect_vector(typify(sq_dna, "dna_ext"),
                ptype = sq_ptype(CPP_get_standard_alphabet("dna_ext"), "dna_ext"),
                size = vec_size(sq_dna))
})

# ERROR FOR NON-SQ OBJECTS ----
test_sq_only(typify, "dna_bsc")

# NO CHANGES IF ALREADY TARGET CLASS ----
test_that("typify() returns unaltered sq when sq is already of target class", {
  expect_identical(typify(sq_dna, "dna_bsc"), sq_dna)
})

# FAIL WHEN OUT-OF-ALPHABET LETTERS ----
test_that("typify() throws error when there are letters not in the target alphabet", {
  expect_error(typify(sq_unt, "dna_bsc"))
})

test_that("typify() treats multiple-character letters as a whole (and throws error)", {
  # If each character was considered separate, character-only sequences like sq_atp
  # would be typified to AMI_EXT, as AMI_EXT alphabet contains all LETTERS
  expect_error(typify(sq_atp, "ami_ext"))
})

test_that("typify() returns object that is convertible to identical character vector ", {
  expect_equal(as.character(typify(sq_dna, "dna_ext")),
               as.character(sq_dna))
  expect_equal(as.character(typify(sq_unt, "ami_ext")),
               as.character(sq_unt))
})

# EDGE CASES ----
test_that("typify() works for empty sq and returns empty sq of target class", {
  expect_vector(typify(sq_empty, "ami_bsc"),
                ptype = sq_ptype(CPP_get_standard_alphabet("ami_bsc"), "ami_bsc"),
                size = vec_size(sq_empty))
})
