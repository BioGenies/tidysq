# SETUP ----
str_dna <- c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT")
sq_dna <- sq(str_dna, "dna_bsc")

fun_1 <- function(sequence) sequence[1]
fun_2 <- function(sequence) sum(sequence == "A")

# ERROR FOR NON-SQ OBJECTS ----
test_that("sqapply() throws an error whenever passed object of class other that sq", {
  expect_error(sqapply(1:7, fun_1))
  expect_error(sqapply(LETTERS, fun_2))
  expect_error(sqapply(list(mean, sum, sd), fun_1))
})

# TYPE CORRECTNESS
test_that("sqapply() returns list", {
  expect_list(sqapply(sq_dna, fun_1))
  expect_list(sqapply(sq_dna, fun_2))
})

test_that("sqapply() with identity function returns list of splitted sequences", {
  expect_identical(sqapply(sq_dna, identity),
                   strsplit(str_dna, ""))
})

test_that("sqapply() with identity function returns list of sequences as strings when single_string = TRUE", {
  expect_identical(sqapply(sq_dna, identity, single_string = TRUE),
                   as.list(str_dna))
})

# APPLICATION ----
test_that("sqapply() returns correct result", {
  expect_equal(sqapply(sq_dna, fun_1),
               list("T", "C", "T", NA_character_, "A"))
  expect_equal(sqapply(sq_dna, fun_2),
               list(2, 2, 2, 0, 1))
})