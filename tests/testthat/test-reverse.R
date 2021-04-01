# SETUP ----
str_dna <- c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT")
str_dna_rev <- c("GTACGGGTCAT", "AGGCTGGAC", "GCCTGATGAT", "", "TGGCA")
str_ami <- c("PADINI", "HOUDINI", "UNBELIEVABLE")
str_ami_rev <- c("INIDAP", "INIDUOH", "ELBAVEILEBNU")

sq_dna <- sq(str_dna, alphabet = "dna_bsc")
sq_ami <- sq(str_ami, alphabet = "ami_ext")

# PROTOTYPE PRESERVATION ----
test_that("reverse() preserves all attributes of original vector", {
  expect_vector(reverse(sq_dna),
                ptype = vec_ptype(sq_dna),
                size = vec_size(sq_dna))
  expect_vector(reverse(sq_ami),
                ptype = vec_ptype(sq_ami),
                size = vec_size(sq_ami))
})

# ERROR FOR NON-SQ OBJECTS ----
test_that("reverse() throws an error whenever passed object of class other that sq", {
  expect_error(reverse(1:7))
  expect_error(reverse(LETTERS))
  expect_error(reverse(list(mean, sum, sd)))
})

# VALUE COMPUTATION ----
test_that("reverse() returns correct value", {
  expect_equal(as.character(reverse(sq_dna)),
               str_dna_rev)
  expect_equal(as.character(reverse(sq_ami)),
               str_ami_rev)
})

# CANCELLING UPON DOUBLE USAGE ----
test_that("double use of reverse() returns original value", {
  expect_identical(reverse(reverse(sq_dna)),
                   sq_dna)
  expect_identical(reverse(reverse(sq_ami)),
                   sq_ami)
})
