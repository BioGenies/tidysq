# SETUP ----
sq_unt <- sq(c("PQNVIXFD", "PDOQXN-FI", "SPBI-F-XXS"), alphabet = "unt")
sq_atp <- sq(c("mAmYmY", "nbAnsAmA", ""),
             alphabet = c("mA", "mY", "nbA", "nsA"))
sq_dna <- sq(c("TTCAGGGCTAG", "CGATTGC", "CAGTTTA"), alphabet = "dna_bsc")
sq_empty <- sq(character(), alphabet = "rna_bsc")

# CORRECT PROTOTYPE OF RETURNED VALUE ----
test_that("typify() returns an sq object of desired type", {
  skip(".apply_sq() inside typify() not implemented")
  expect_vector(typify(sq_unt, "ami_ext"),
                ptype = sq_ptype(get_standard_alphabet("ami_ext"), "ami_ext"),
                size = vec_size(sq_unt))
})

# NO CHANGES IF ALREADY TARGET CLASS ----
test_that("typify() returns unaltered sq when sq is already of target class", {
  skip(".apply_sq() inside typify() not implemented")
  expect_reference(typify(sq_dna, "dna_bsc"), sq_dna)
})

# FAIL WHEN OUT-OF-ALPHABET LETTERS ----
test_that("typify() throws error when there are letters not in the target alphabet", {
  skip(".apply_sq() inside typify() not implemented")
  expect_error(typify(sq_unt, "dna_bsc"))
})

test_that("typify() treats multiple-character letters as a whole (and throws error)", {
  # If each character was considered separate, character-only sequences like sq_atp
  # would be typified to amisq, as amisq alphabet contains all LETTERS
  skip(".apply_sq() inside typify() not implemented")
  expect_error(typify(sq_atp, "ami_ext"))
})

# EDGE CASES ----
test_that("typify() works for empty sq and returns empty sq of target class", {
  skip(".apply_sq() inside typify() not implemented")
  expect_vector(typify(sq_empty, "ami_bsc"),
                ptype = sq_ptype(get_standard_alphabet("ami_ext"), "ami_bsc"),
                size = vec_size(sq_empty))
})
