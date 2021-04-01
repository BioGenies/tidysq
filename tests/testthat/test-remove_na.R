# SETUP ----
suppressWarnings({
  sq_ami <- bite(sq(c("MIAANYTWIL", "", "TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"),
                    alphabet = "ami_ext"), 1:14)
})
sq_ami_2 <- sq(c("MIAANYTWIL", "", "TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"),
               alphabet = "ami_ext")
sq_ami_3 <- sq(c("", "", "TIAALGNIIYRAIE", "", ""),
               alphabet = "ami_ext")
suppressWarnings({
  sq_dna <- bite(sq(c("ATGCAGGA", "GACCGAACGAN", "TGACGAGCTTA", "ACTNNAGCN"),
                    alphabet = "dna_ext"), 1:11)
})
sq_dna_2 <- sq(c("ATGCAGGA", "GACCGAACGAN", "TGACGAGCTTA", "ACTNNAGCN"),
               alphabet = "dna_ext")
sq_dna_3 <- sq(c("", "GACCGAACGAN", "TGACGAGCTTA", ""),
               alphabet = "dna_ext")
suppressWarnings({
  sq_rna <- bite(sq(c("UAAACGGGCUA", "UACACGGGCGACU", "AGGCA", "AUGGCGGAUGUUC"),
                    alphabet = "rna_bsc"), 1:13)
})
sq_rna_2 <- sq(c("UAAACGGGCUA", "UACACGGGCGACU", "AGGCA", "AUGGCGGAUGUUC"),
               alphabet = "rna_bsc")
sq_rna_3 <- sq(c("", "UACACGGGCGACU", "", "AUGGCGGAUGUUC"),
               alphabet = "rna_bsc")

# PROTOTYPE PRESERVATION ----
test_that("remove_na() with `by_letter = TRUE` preserves all attributes of original vector", {
  expect_vector(remove_na(sq_ami, by_letter = TRUE),
                ptype = vec_ptype(sq_ami),
                size = vec_size(sq_ami))
  expect_vector(remove_na(sq_dna, by_letter = TRUE),
                ptype = vec_ptype(sq_dna),
                size = vec_size(sq_dna))
  expect_vector(remove_na(sq_rna, by_letter = TRUE),
                ptype = vec_ptype(sq_rna),
                size = vec_size(sq_rna))
})

test_that("remove_na() with `by_letter = FALSE` preserves all attributes of original vector", {
  expect_vector(remove_na(sq_ami, by_letter = FALSE),
                ptype = vec_ptype(sq_ami),
                size = vec_size(sq_ami))
  expect_vector(remove_na(sq_dna, by_letter = FALSE),
                ptype = vec_ptype(sq_dna),
                size = vec_size(sq_dna))
  expect_vector(remove_na(sq_rna, by_letter = FALSE),
                ptype = vec_ptype(sq_rna),
                size = vec_size(sq_rna))
})

# ERROR FOR NON-SQ OBJECTS ----
test_that("remove_na() throws an error whenever passed object of class other that sq", {
  expect_error(remove_na(1:7))
  expect_error(remove_na(LETTERS))
  expect_error(remove_na(list(mean, sum, sd)))
})

# VALUE COMPUTATION ----
test_that("remove_na() with `by_letter = TRUE` returns correct value", {
  expect_equal(as.character(remove_na(sq_ami, by_letter = TRUE)),
               as.character(sq_ami_2))
  expect_equal(as.character(remove_na(sq_dna, by_letter = TRUE)),
               as.character(sq_dna_2))
  expect_equal(as.character(remove_na(sq_rna, by_letter = TRUE)),
               as.character(sq_rna_2))
})

test_that("remove_na() with `by_letter = FALSE` returns correct value", {
  expect_equal(as.character(remove_na(sq_ami, by_letter = FALSE)),
               as.character(sq_ami_3))
  expect_equal(as.character(remove_na(sq_dna, by_letter = FALSE)),
               as.character(sq_dna_3))
  expect_equal(as.character(remove_na(sq_rna, by_letter = FALSE)),
               as.character(sq_rna_3))
})

# NO EFFECT ON SEQUENCES WITHOUT NA'S ----
test_that("remove_na() has no effect when original value has no NA's", {
  expect_identical(remove_na(sq_ami_2, by_letter = TRUE),
                   sq_ami_2)
  expect_identical(remove_na(sq_ami_2, by_letter = FALSE),
                   sq_ami_2)
})
