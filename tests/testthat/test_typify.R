# SETUP ----
sq_unt <- .construct_untsq(c("PQNVIXFD", "PDOQXN-FI", "SPBI-F-XXS"))
sq_atp <- .nonst_construct_sq(c("mAmYmY", "nbAnsAmA", ""), c("mA", "mY", "nbA", "nsA"))
sq_dna <- construct_sq_dna(c("TTCAGGGCTAG", "CGATTGC", "CAGTTTA"), is_clean = TRUE)
sq_empty <- construct_sq_rna(character(), is_clean = TRUE)

# CORRECT PROTOTYPE OF RETURNED VALUE ----
test_that("typify() returns an sq object of desired type", {
  expect_vector(typify(sq_unt, "ami"),
                ptype = .construct_sq_ptype("ami", FALSE),
                size = vec_size(sq_unt))
})

# NO CHANGES IF ALREADY TARGET CLASS ----
test_that("typify() returns unaltered sq when sq is already of target class", {
  expect_reference(typify(sq_dna, "dna"), sq_dna)
})

# FAIL WHEN OUT-OF-ALPHABET LETTERS ----
test_that("typify() throws error when there are letters not in the target alphabet", {
  expect_error(typify(sq_unt, "dna"))
})

test_that("typify() treats multiple-character letters as a whole (and throws error)", {
  # If each character was considered separate, character-only sequences like sq_atp
  # would be typified to amisq, as amisq alphabet contains all LETTERS
  expect_error(typify(sq_atp, "ami"))
})

# EDGE CASES ----
test_that("typify() works for empty sq and returns empty sq of target class", {
  expect_vector(typify(sq_empty, "ami"),
                ptype = .construct_sq_ptype("ami", FALSE),
                size = vec_size(sq_empty))
})
