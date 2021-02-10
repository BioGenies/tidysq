# SETUP ----
sq_ami <- sq(c("MIAANYTWIL","TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"),
             alphabet = "ami_ext")
sq_ami_2 <- sq(c("MIAANYTWIL","TIAALGNIIYRAIE", "NYERTGHLI", ""),
               alphabet = "ami_bsc")
sq_ami_3 <- sq(c("MIAANYTWIL","TIAALGNIIYRAIE", "NYERTGHLI", "MAYIALN"),
               alphabet = "ami_bsc")
sq_dna <- sq(c("ATGCAGGA", "GACCGAACGAN", "TGACGAGCTTA", "ACTNNAGCN"),
             alphabet = "dna_ext")
sq_dna_2 <- sq(c("ATGCAGGA", "", "TGACGAGCTTA", ""),
               alphabet = "dna_bsc")
sq_dna_3 <- sq(c("ATGCAGGA", "GACCGAACGA", "TGACGAGCTTA", "ACTAGC"),
               alphabet = "dna_bsc")
sq_rna <- sq(c("UGGCGNNGBV", "ACGGUUUCGUU", "UBBGDGAACG", "GGCUCGACAGACUGC"),
             alphabet = "rna_ext")
sq_rna_2 <- sq(c("", "ACGGUUUCGUU", "", "GGCUCGACAGACUGC"),
               alphabet = "rna_bsc")
sq_rna_3 <- sq(c("UGGCGG", "ACGGUUUCGUU", "UGGAACG", "GGCUCGACAGACUGC"),
               alphabet = "rna_bsc")

# CORRECT PROTOTYPE OF RETURNED VALUE ----
test_that("remove_ambiguous() returns an sq object with _bsc class", {
  expect_vector(remove_ambiguous(sq_ami),
                ptype = sq_ptype(CPP_get_standard_alphabet("ami_bsc"), "ami_bsc"),
                size = vec_size(sq_ami))
  expect_vector(remove_ambiguous(sq_dna_3),
                ptype = sq_ptype(CPP_get_standard_alphabet("dna_bsc"), "dna_bsc"),
                size = vec_size(sq_dna_3))
  expect_vector(remove_ambiguous(sq_rna_2),
                ptype = sq_ptype(CPP_get_standard_alphabet("rna_bsc"), "rna_bsc"),
                size = vec_size(sq_rna_2))
})

# ERROR FOR NON-STANDARD SQ OBJECTS ----
test_that("remove_ambiguous() throws an error whenever passed object of class other that standard sq classes", {
  expect_error(remove_ambiguous(1:7))
  expect_error(remove_ambiguous(LETTERS))
  expect_error(remove_ambiguous(list(mean, sum, sd)))
  expect_error(remove_ambiguous(sq(c(")R#)#!Vawr9fy", "*V)RUgBa^%#!b]"))))
  expect_error(remove_ambiguous(sq(c("accmsce", "auprcacc"), alphabet = c("auprc", "acc", "msce"))))
})

# VALUE COMPUTATION ----
test_that("remove_ambiguous() removes whole sequences correctly", {
  expect_equal(remove_ambiguous(sq_ami),
               sq_ami_2)
  expect_equal(remove_ambiguous(sq_dna),
               sq_dna_2)
  expect_equal(remove_ambiguous(sq_rna),
               sq_rna_2)
})

test_that("remove_ambiguous() removes whole sequences when `by_letter = FALSE`", {
  expect_equal(remove_ambiguous(sq_ami, by_letter = FALSE),
               sq_ami_2)
  expect_equal(remove_ambiguous(sq_dna, by_letter = FALSE),
               sq_dna_2)
  expect_equal(remove_ambiguous(sq_rna, by_letter = FALSE),
               sq_rna_2)
})

test_that("remove_ambiguous() removes letters only when `by_letter = TRUE`", {
  expect_equal(remove_ambiguous(sq_ami, by_letter = TRUE),
               sq_ami_3)
  expect_equal(remove_ambiguous(sq_dna, by_letter = TRUE),
               sq_dna_3)
  expect_equal(remove_ambiguous(sq_rna, by_letter = TRUE),
               sq_rna_3)
})

# NO CHANGES IF ALREADY TARGET CLASS ----
test_that("remove_ambiguous() doesn't affect basic alphabet sequences", {
  expect_identical(remove_ambiguous(sq_ami_2), sq_ami_2)
  expect_identical(remove_ambiguous(sq_dna_3), sq_dna_3)
  expect_identical(remove_ambiguous(sq_rna_2), sq_rna_2)
})
