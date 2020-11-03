# SETUP ----
str_dna <- c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT")
str_rna <- c("", "KBS-UVW-AWWWG", "YGHHH-", "-CRASH", "MND-KUUBV-MY-")
str_ami <- c("OUTLANDISH", "UNSTRUCTURIZED", "FEAR")
str_unt <- c("vip01", "vip02", "vip04", "missing_one")
str_atp <- c("mAmYmY", "nbAnsAmA", "")

dna_bsc_alph <- c("A", "C", "G", "T", "-")
rna_ext_alph <- c("A", "C", "G", "U", "W", "S", "M", "K", "R", "Y", "B", "D", 
                  "H", "V", "N", "-")

ami_ext_alph <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", 
                  "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", 
                  "Z", "-", "*")
atp_alph <- c("mA", "mY", "nbA", "nsA")

# CORRECT PROTOTYPE OF RETURNED VALUE ----
test_that("construct_sq() returns object of correct prototype", {
  expect_vector(sq(str_dna, alphabet = "dna_bsc"),
                ptype = sq_ptype(dna_bsc_alph, "dna_bsc"),
                size = vec_size(str_dna))
  expect_vector(sq(str_rna, alphabet = "rna_ext"),
                ptype = sq_ptype(rna_ext_alph, "rna_ext"),
                size = vec_size(str_rna))
  expect_vector(sq(str_ami, alphabet = "ami_ext"),
                ptype = sq_ptype(ami_ext_alph, "ami_ext"),
                size = vec_size(str_ami))
})

test_that("construct_sq() returns object of correct class for unt and atp options", {
  expect_s3_class(sq(str_atp, alphabet = atp_alph),
                  class = "sq_atp",
                  exact = FALSE)
  
  skip("untyped not implemented yet")
  expect_s3_class(sq(str_unt, alphabet = "unt"),
                  class = "sq_unt",
                  exact = FALSE)
})

test_that("construct_sq() returns object of same size as passed character vector", {
  expect_equal(vec_size(sq(str_atp, alphabet = atp_alph)),
               vec_size(str_atp))
  
  skip("untyped not implemented yet")
  expect_equal(vec_size(sq(str_unt, alphabet = "unt")),
               vec_size(str_unt))
})

# TODO: test C_get_real_alph() somewhere
test_that("construct_sq() returns object with alphabet attribute that contains existing letters for unt and atp options", {
  expect_setequal(
    alphabet(sq(str_atp, alphabet = atp_alph)),
    atp_alph
  )
  
  skip("untyped not implemented yet")
  expect_setequal(
    alphabet(sq(str_unt, alphabet = "unt")),
    C_get_real_alph(str_unt)
  )
})

# REVERSING TO CHARACTER ----
test_that("applying to.character() returns original character vector", {
  expect_equivalent(
    as.character(sq(str_dna, alphabet = "dna_bsc")),
    str_dna
  )
  expect_equivalent(
    as.character(sq(str_rna, alphabet = "rna_ext")),
    str_rna
  )
  expect_equivalent(
    as.character(sq(str_ami, alphabet = "ami_ext")),
    str_ami
  )
  expect_equivalent(
    as.character(sq(str_atp, alphabet = atp_alph)),
    str_atp
  )
  
  skip("untyped not implemented yet")
  expect_equivalent(
    as.character(sq(str_unt, alphabet = "unt")),
    str_unt
  )
})


# SHORTHAND FUNCTIONS ----
# TODO: implement type guessing 
test_that("construct_sq_dna() return identical value as construct_sq() with \"dna\" type", {
  skip("type guessing not implemented yet")
  expect_identical(construct_sq_dna(str_dna, is_clean = TRUE),
                   construct_sq(str_dna, type = "dna", is_clean = TRUE))
})
test_that("construct_sq_rna() return identical value as construct_sq() with \"rna\" type", {
  skip("type guessing not implemented yet")
  expect_identical(construct_sq_rna(str_rna, is_clean = FALSE),
                   construct_sq(str_rna, type = "rna", is_clean = FALSE))
})
test_that("construct_sq_ami() return identical value as construct_sq() with \"ami\" type", {
  skip("type guessing not implemented yet")
  expect_identical(construct_sq_ami(str_ami, is_clean = FALSE),
                   construct_sq(str_ami, type = "ami", is_clean = FALSE))
})
