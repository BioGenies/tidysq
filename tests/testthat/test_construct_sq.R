# SETUP ----
str_dna <- c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT")
str_rna <- c("", "KBS-UVW-AWWWG", "YGHHH-", "-CRASH", "MND-KUUBV-MY-")
str_ami <- c("OUTLANDISH", "UNSTRUCTURIZED", "FEAR")
str_unt <- c("vip01", "vip02", "vip04", "missing_one")
str_atp <- c("mAmYmY", "nbAnsAmA", "")

# CORRECT PROTOTYPE OF RETURNED VALUE ----
test_that("construct_sq() returns object of correct prototype", {
  expect_vector(construct_sq(str_dna, type = "dna", is_clean = TRUE),
                ptype = .construct_sq_ptype("dna", is_clean = TRUE),
                size = vec_size(str_dna))
  expect_vector(construct_sq(str_rna, type = "rna", is_clean = FALSE),
                ptype = .construct_sq_ptype("rna", is_clean = FALSE),
                size = vec_size(str_rna))
  expect_vector(construct_sq(str_ami, type = "ami", is_clean = FALSE),
                ptype = .construct_sq_ptype("ami", is_clean = FALSE),
                size = vec_size(str_ami))
})

test_that("construct_sq() returns object of correct class for unt and atp options", {
  expect_s3_class(construct_sq(str_unt, type = "unt"),
                  class = "untsq",
                  exact = FALSE)
  expect_s3_class(construct_sq(str_atp, non_standard = c("mA", "mY", "nbA", "nsA")),
                  class = "atpsq",
                  exact = FALSE)
})

test_that("construct_sq() returns object of same size as passed character vector", {
  expect_equal(vec_size(construct_sq(str_unt, type = "unt")),
               vec_size(str_unt))
  expect_equal(vec_size(construct_sq(str_atp, non_standard = c("mA", "mY", "nbA", "nsA"))),
               vec_size(str_atp))
})

# TODO: test C_get_real_alph() somewhere
test_that("construct_sq() returns object with alphabet attribute that contains existing letters for unt and atp options", {
  expect_setequal(
    alphabet(construct_sq(str_unt, type = "unt")),
    C_get_real_alph(str_unt)
  )
  expect_setequal(
    alphabet(construct_sq(str_atp, non_standard = c("mA", "mY", "nbA", "nsA"))),
    c("mA", "mY", "nbA", "nsA")
  )
})

# REVERSING TO CHARACTER ----
test_that("applying to.character() returns original character vector", {
  expect_equivalent(
    as.character(construct_sq(str_dna, type = "dna", is_clean = TRUE)),
    str_dna
  )
  expect_equivalent(
    as.character(construct_sq(str_rna, type = "rna", is_clean = FALSE)),
    str_rna
  )
  expect_equivalent(
    as.character(construct_sq(str_ami, type = "ami", is_clean = FALSE)),
    str_ami
  )
  expect_equivalent(
    as.character(construct_sq(str_unt, type = "unt")),
    str_unt
  )
  expect_equivalent(
    as.character(construct_sq(str_atp, non_standard = c("mA", "mY", "nbA", "nsA"))),
    str_atp
  )
})

# SHORTHAND FUNCTIONS ----
test_that("construct_sq_dna() return identical value as construct_sq() with \"dna\" type", {
  expect_identical(construct_sq_dna(str_dna, is_clean = TRUE),
                   construct_sq(str_dna, type = "dna", is_clean = TRUE))
})
test_that("construct_sq_rna() return identical value as construct_sq() with \"rna\" type", {
  expect_identical(construct_sq_rna(str_rna, is_clean = FALSE),
                   construct_sq(str_rna, type = "rna", is_clean = FALSE))
})
test_that("construct_sq_ami() return identical value as construct_sq() with \"ami\" type", {
  expect_identical(construct_sq_ami(str_ami, is_clean = FALSE),
                   construct_sq(str_ami, type = "ami", is_clean = FALSE))
})
