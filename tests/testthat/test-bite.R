# SETUP ----
char_atp <- c("mA", "mY", "nbA", "nsA")

sq_unt <- sq(c("PQNVIXFD", "PDOQXN-FI", "SPBI-F-XXS"), alphabet = "unt")
sq_atp <- sq(c("mAmYmY", "nbAnsAmA", ""), alphabet = char_atp)
sq_dna <- sq(c("TTCAGGGCTAG", "CGATTGC", "CAGTTTA"), alphabet = "dna_bsc")
sq_empty <- sq(character(), alphabet = "rna_bsc")

NA_warning <- "some sequences are subsetted with index bigger than length - NA introduced"
NA_letter <- getOption("tidysq_NA_letter")

# CORRECT PROTOTYPE OF RETURNED VALUE ----
test_that("bite() returns an sq object of original type", {
  suppressWarnings({
    expect_vector(bite(sq_unt, 2:5),
                  ptype = sq_ptype(vec_data(alphabet(sq_unt)), "unt"),
                  size = vec_size(sq_unt))
    expect_vector(bite(sq_atp, 2:5),
                  ptype = sq_ptype(char_atp, "atp"),
                  size = vec_size(sq_atp))
    expect_vector(bite(sq_dna, 2:5),
                  ptype = sq_ptype(CPP_get_standard_alphabet("dna_bsc"), "dna_bsc"),
                  size = vec_size(sq_dna))
    expect_vector(bite(sq_empty, 2:5),
                  ptype = sq_ptype(CPP_get_standard_alphabet("rna_bsc"), "rna_bsc"),
                  size = vec_size(sq_empty))
  })
})

# ERROR FOR NON-SQ OBJECTS ----
test_that("bite() throws an error whenever passed object of class other that sq", {
  expect_error(bite(1:7, 4:6))
  expect_error(bite(LETTERS, -7:-11))
  expect_error(bite(list(mean, sum, sd), 1))
})

# HANDLING INDICES INSIDE SEQUENCE ORIGINAL LENGTH ----
test_that("bite() interprets positive indices correctly", {
  expect_equivalent(as.character(bite(sq_unt, c(4, 2, 1))),
                    c("VQP", "QDP", "IPS"))
  expect_equivalent(as.character(bite(sq_dna, 4:7)),
                    c("AGGG", "TTGC", "TTTA"))
})

test_that("bite() correctly interpretes '-1' index (unknown as '-0' in C++)", {
  expect_equivalent(as.character(bite(sq_dna, -1)),
                    c("TCAGGGCTAG", "GATTGC", "AGTTTA"))
})

test_that("bite() interprets negative indices correctly", {
  expect_equivalent(as.character(bite(sq_unt, -6)),
                    c("PQNVIFD", "PDOQX-FI", "SPBI--XXS"))
  expect_equivalent(as.character(bite(sq_dna, -3:-7)),
                    c("TTCTAG", "CG", "CA"))
})

# HANDLING INDICES OUT OF SEQUENCE ORIGINAL LENGTH ----
test_that("bite() raises a warning when indices reach outside original length", {
  expect_warning(bite(sq_unt, 13:16), NA_warning)
  expect_warning(bite(sq_atp, 13:16), NA_warning)
  expect_warning(bite(sq_dna, 13:16), NA_warning)
})

test_that("bite() returns NA_letter for each index outside of original length", {
  suppressWarnings({
    expect_equivalent(as.character(bite(sq_unt, 13:16)),
                      rep(paste0(rep(NA_letter, 4), collapse = ""), vec_size(sq_unt)))
    expect_equivalent(as.character(bite(sq_atp, 3:5)),
                      paste0(c("mY", "mA", NA_letter),
                             paste0(rep(NA_letter, 2), collapse = "")))
    expect_equivalent(as.character(bite(sq_dna, c(16, 20, 13))),
                      rep(paste0(rep(NA_letter, 3), collapse = ""), vec_size(sq_dna)))
  })
})

test_that("bite() ignores negative indices outside of original length", {
  expect_equivalent(as.character(bite(sq_unt, c(-6, -20))),
                    c("PQNVIFD", "PDOQX-FI", "SPBI--XXS"))
  expect_equivalent(as.character(bite(sq_dna, -15:-19)),
                    c("TTCAGGGCTAG", "CGATTGC", "CAGTTTA"))
  expect_equivalent(as.character(bite(sq_atp, c(-2, -5))),
                    c("mAmY", "nbAmA", ""))
})

# CORNER/EDGE CASES ----
test_that("bite() does nothing on empty sq objects", {
  expect_equal(bite(sq_empty, 2), sq_empty)
  expect_equal(bite(sq_empty, 5:9), sq_empty)
  expect_equal(bite(sq_empty, c(7, 16, 3)), sq_empty)
  expect_equal(bite(sq_empty, -c(8, 2, 1)), sq_empty)
})

test_that("bite() ignores multiple instances of the same negative index", {
  expect_identical(as.character(bite(sq_unt, c(-6, -6, -6))),
                   c("PQNVIFD", "PDOQX-FI", "SPBI--XXS"))
  expect_identical(as.character(bite(sq_dna, c(-4, -1, -1, -7, -4))),
                   c("TCGGCTAG", "GATG", "AGTT"))
})

test_that("bite() throws an error when passed mixed positive and negative indices", {
  expect_error(bite(sq_dna, -5:5))
  expect_error(bite(sq_empty, c(6, -1, 0, 5)))
})
