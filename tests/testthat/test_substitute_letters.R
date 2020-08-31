# SETUP ----
sq_dna <- construct_sq_dna(c("ACTGTC", "CGCGTTA"), is_clean = TRUE)
sq_ami <- construct_sq_ami(c("APOPNIQEV", "CSVMIBF"), is_clean = FALSE)
# sq_unt <- construct_sq(c("GO%NC@E(123)RO", "NFI%(#)VT;"))

# CORRECT PROTOTYPE OF RETURNED VALUE ----
test_that("substitute_letters() returns an atpsq object", {
  expect_s3_class(substitute_letters(sq_dna, c(T = "U")),
                  class = "atpsq",
                  exact = FALSE)
  expect_s3_class(substitute_letters(sq_ami, c(P = "Q", O = "W", I = "E", U = "R", Y = "T", G = "V",
                                               L = "A", K = "S", J = "D", H = "F", M = "Z", N = "X", B = "C")),
                  class = "atpsq",
                  exact = FALSE)
})

test_that("substitute_letters() returns an object with trimmed alphabet attribute", {
  expect_setequal(
    alphabet(substitute_letters(sq_dna, c(T = "U"))),
    c("A", "C", "G", "U", "-")
  )
  expect_setequal(
    alphabet(substitute_letters(sq_ami, c(P = "Q", O = "W", I = "E", U = "R", Y = "T", G = "V",
                                          L = "A", K = "S", J = "D", H = "F", M = "Z", N = "X", B = "C"))),
    .skip_characters(alphabet(sq_ami), c("P", "O", "I", "U", "Y", "G", "L", "K", "J", "H", "M", "N", "B"))
  )
})

test_that("substitute_letters() does not remove nor add sequences", {
  expect_length(substitute_letters(sq_dna, c(T = "U")),
                length(sq_dna))
  expect_length(substitute_letters(sq_ami, c(P = "Q", O = "W", I = "E", U = "R", Y = "T", G = "V",
                                             L = "A", K = "S", J = "D", H = "F", M = "Z", N = "X", B = "C")),
                length(sq_ami))
})

test_that("substitute_letters() keep original_lengths unchanged", {
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
