# SETUP ----
sq_dna <- construct_sq_dna(c("TTCAGGGCTAG", "CGATTGC", "CAGTTTA"), is_clean = TRUE)
sq_rna <- construct_sq_rna(c("UAGUAACCGUAAGCG", "UAGUCC--UA-G"), is_clean = TRUE)

# CORRECT PROTOTYPE OF RETURNED VALUE ----
test_that("encode() returns an encsq object", {
  expect_s3_class(encode(sq_dna, c(A = 1, C = 2, G = 3, T = 4, `-` = 7)),
                  class = "encsq",
                  exact = FALSE)
  expect_s3_class(encode(sq_rna, c(A = 1, C = 2, G = 3, U = 20)),
                  class = "encsq",
                  exact = FALSE)
})

test_that("encode() returns value with numeric alphabet attribute", {
  expect_setequal(
    alphabet(encode(sq_dna, c(A = 1, C = 2, G = 3, T = 4, `-` = 7))),
    c(1, 2, 3, 4, 7)
  )
})

test_that("numeric alphabet has NAs whenever no encoding is given for a letter", {
  expect_setequal(
    alphabet(encode(sq_rna, c(A = 1, C = 2, G = 3, U = 20))),
    c(1, 2, 3, 20, NA_real_)
  )
})

test_that("encode() does not remove nor add sequences", {
  expect_length(encode(sq_dna, c(A = 1, C = 2, G = 3, T = 4, `-` = 7)),
                length(sq_dna))
  expect_length(encode(sq_rna, c(A = 1, C = 2, G = 3, U = 20)),
                length(sq_rna))
})

test_that("encode() keep original_lengths unchanged", {
  expect_equal(
    get_sq_lengths(encode(sq_dna, c(A = 1, C = 2, G = 3, T = 4, `-` = 7))),
    get_sq_lengths(sq_dna)
  )
  expect_equal(
    get_sq_lengths(encode(sq_rna, c(A = 1, C = 2, G = 3, U = 20))),
    get_sq_lengths(sq_rna)
  )
})

# PRESERVATION OF INTERNAL REPRESENTATION ----
test_that("encode() only alters attributes of sq object", {
  expect_identical(
    vec_data(encode(sq_dna, c(A = 1, C = 2, G = 3, T = 4, `-` = 7))),
    vec_data(sq_dna)
  )
  expect_identical(
    vec_data(encode(sq_rna, c(A = 1, C = 2, G = 3, U = 20))),
    vec_data(sq_rna)
  )
})
