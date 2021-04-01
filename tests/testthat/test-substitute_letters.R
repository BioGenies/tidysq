# SETUP ----
sq_dna <- sq(c("ACTGTC", "CGCGTTA"), alphabet = "dna_bsc")
sq_ami <- sq(c("APOPNIQEV", "CSVMIBF"), alphabet = "ami_ext")
# sq_unt <- sq(c("GO%NC@E(123)RO", "NFI%(#)VT;"), alphabet = "unt")

str_dna <- c("ACBGBC", "CGCGBBA")
str_dna_num <- c("AC6G6C", "CGCG66A")
str_ami <- c("AQWQXEQEV", "CSVZECF")
str_ami_num <- c("A11O11NIQ3322", "CS22MI7F")

# CORRECT PROTOTYPE OF RETURNED VALUE ----
test_that("substitute_letters() returns an sq_atp object", {
  expect_s3_class(substitute_letters(sq_dna, c(T = "B")),
                  class = "sq_atp",
                  exact = FALSE)
  expect_s3_class(substitute_letters(sq_ami, c(P = "Q", O = "W", I = "E", U = "R", Y = "T", G = "V",
                                               L = "A", K = "S", J = "D", H = "F", M = "Z", N = "X", B = "C")),
                  class = "sq_atp",
                  exact = FALSE)
})

test_that("substitute_letters() returns an object with trimmed alphabet attribute", {
  expect_setequal(
    alphabet(substitute_letters(sq_dna, c(T = "B"))),
    c("A", "C", "G", "B", "-")
  )
  expect_setequal(
    alphabet(substitute_letters(sq_ami, c(P = "Q", O = "W", I = "E", U = "R", Y = "T", G = "V",
                                          L = "A", K = "S", J = "D", H = "F", M = "Z", N = "X", B = "C"))),
    setdiff(alphabet(sq_ami), c("P", "O", "I", "U", "Y", "G", "L", "K", "J", "H", "M", "N", "B"))
  )
})

test_that("substitute_letters() does not remove nor add sequences", {
  expect_length(substitute_letters(sq_dna, c(T = "B")),
                length(sq_dna))
  expect_length(substitute_letters(sq_ami, c(P = "Q", O = "W", I = "E", U = "R", Y = "T", G = "V",
                                             L = "A", K = "S", J = "D", H = "F", M = "Z", N = "X", B = "C")),
                length(sq_ami))
})

test_that("substitute_letters() keep original_lengths unchanged", {
  expect_equal(
    get_sq_lengths(substitute_letters(sq_dna, c(T = "B"))),
    get_sq_lengths(sq_dna)
  )
  expect_equal(
    get_sq_lengths(substitute_letters(sq_ami, c(P = "Q", O = "W", I = "E", U = "R", Y = "T", G = "V",
                                                L = "A", K = "S", J = "D", H = "F", M = "Z", N = "X", B = "C"))),
    get_sq_lengths(sq_ami)
  )
})

# ARGUMENT PREREQUISITES ----
test_that("substitute_letters() throws an error whenever passed object of class other that sq", {
  expect_error(substitute_letters(1:7, c(S = "H")))
  expect_error(substitute_letters(LETTERS, c(S = "SH", H = "HS")))
  expect_error(substitute_letters(list(mean, sum, sd), c(`mean` = "Aix")))
})

test_that("encoding must be either numeric or character", {
  expect_error(substitute_letters(sq_dna, c(A = TRUE, T = FALSE)),
               "encoding must be either numeric of character vector")
  expect_error(substitute_letters(sq_dna, list(G = mean, C = sum)))
})

# CORRECT RETURN VALUE ----
test_that("substitute_letters() correctly computes value", {
  expect_equal(
    as.character(substitute_letters(sq_dna, c(T = "B"))),
    str_dna
  )
  expect_equal(
    as.character(substitute_letters(sq_ami, c(P = "Q", O = "W", I = "E", U = "R", Y = "T", G = "V",
                                              L = "A", K = "S", J = "D", H = "F", M = "Z", N = "X", B = "C"))),
    str_ami
  )
})

# COERCING NUMERIC VECTORS TO CHARACTERS ----
test_that("substitute_letters() coerces numeric encodings to string vectors", {
  expect_equal(
    as.character(substitute_letters(sq_dna, c(T = 6))),
    str_dna_num
  )
  expect_equal(
    as.character(substitute_letters(sq_ami, c(P = 11, V = 22, E = 33, X = 6, B = 7))),
    str_ami_num
  )
})

# SURJECTIVITY OF SUBSTITUTION ----
test_that("substitute_letters() is a surjection regarding the alphabets", {
  expect_lte(
    vec_size(alphabet(substitute_letters(sq_dna, c(T = "B")))),
    vec_size(alphabet(sq_dna))
  )
  expect_lte(
    vec_size(alphabet(substitute_letters(sq_ami, c(P = "Q", O = "W", I = "E", U = "R", Y = "T", G = "V",
                                                   L = "A", K = "S", J = "D", H = "F", M = "Z", N = "X", B = "C")))),
    vec_size(alphabet(sq_ami))
  )
})
