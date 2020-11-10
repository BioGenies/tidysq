# SETUP ----
sq_dna <- sq(c("ACTGTC", "CGCGTTA"), alphabet = "dna_bsc")
sq_ami <- sq(c("APOPNIQEV", "CSVMIBF"), alphabet = "ami_ext")
# sq_unt <- sq(c("GO%NC@E(123)RO", "NFI%(#)VT;"), alphabet = "unt")

# CORRECT PROTOTYPE OF RETURNED VALUE ----
test_that("substitute_letters() returns an sq_atp object", {
  skip(".apply_sq() inside substitute_letters() not implemented")
  expect_s3_class(substitute_letters(sq_dna, c(T = "U")),
                  class = "sq_atp",
                  exact = FALSE)
  expect_s3_class(substitute_letters(sq_ami, c(P = "Q", O = "W", I = "E", U = "R", Y = "T", G = "V",
                                               L = "A", K = "S", J = "D", H = "F", M = "Z", N = "X", B = "C")),
                  class = "sq_atp",
                  exact = FALSE)
})

test_that("substitute_letters() returns an object with trimmed alphabet attribute", {
  skip(".apply_sq() inside substitute_letters() not implemented")
  expect_setequal(
    alphabet(substitute_letters(sq_dna, c(T = "U"))),
    c("A", "C", "G", "U", "-")
  )
  expect_setequal(
    alphabet(substitute_letters(sq_ami, c(P = "Q", O = "W", I = "E", U = "R", Y = "T", G = "V",
                                          L = "A", K = "S", J = "D", H = "F", M = "Z", N = "X", B = "C"))),
    setdiff(alphabet(sq_ami), c("P", "O", "I", "U", "Y", "G", "L", "K", "J", "H", "M", "N", "B"))
  )
})

test_that("substitute_letters() does not remove nor add sequences", {
  skip(".apply_sq() inside substitute_letters() not implemented")
  expect_length(substitute_letters(sq_dna, c(T = "U")),
                length(sq_dna))
  expect_length(substitute_letters(sq_ami, c(P = "Q", O = "W", I = "E", U = "R", Y = "T", G = "V",
                                             L = "A", K = "S", J = "D", H = "F", M = "Z", N = "X", B = "C")),
                length(sq_ami))
})

test_that("substitute_letters() keep original_lengths unchanged", {
  skip(".apply_sq() inside substitute_letters() not implemented")
  expect_equal(
    get_sq_lengths(substitute_letters(sq_dna, c(T = "U"))),
    get_sq_lengths(sq_dna)
  )
  expect_equal(
    get_sq_lengths(substitute_letters(sq_ami, c(P = "Q", O = "W", I = "E", U = "R", Y = "T", G = "V",
                                                L = "A", K = "S", J = "D", H = "F", M = "Z", N = "X", B = "C"))),
    get_sq_lengths(sq_ami)
  )
})

# SURJECTIVITY OF SUBSTITUTION ----
test_that("substitute_letters() is a surjection regarding the alphabets", {
  skip(".apply_sq() inside substitute_letters() not implemented")
  expect_lte(
    vec_size(alphabet(substitute_letters(sq_dna, c(T = "U")))),
    vec_size(alphabet(sq_dna))
  )
  expect_lte(
    vec_size(alphabet(substitute_letters(sq_ami, c(P = "Q", O = "W", I = "E", U = "R", Y = "T", G = "V",
                                                   L = "A", K = "S", J = "D", H = "F", M = "Z", N = "X", B = "C")))),
    vec_size(alphabet(sq_ami))
  )
})
