# SETUP ----
sq_ami <- bite(construct_sq(c("MIAANYTWIL", "", "TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"), type = "ami"), 1:14)
sq_ami_no_na_elem <- construct_sq(c("MIAANYTWIL", "", "TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"))
sq_ami_no_na <- construct_sq(c("", "", "TIAALGNIIYRAIE", "", ""))
sq_dna <- bite(construct_sq(c("ATGCAGGA", "GACCGAACGAN", "TGACGAGCTTA", "ACTNNAGCN"), type = "dna"), 1:11)
sq_dna_no_na_elem <- construct_sq(c("ATGCAGGA", "GACCGAACGAN", "TGACGAGCTTA", "ACTNNAGCN"), type = "dna")
sq_dna_no_na <- construct_sq(c("", "GACCGAACGAN", "TGACGAGCTTA", ""), type = "dna")
sq_rna <- bite(construct_sq(c("UAAACGGGCUA", "UACACGGGCGACU", "AGGCA", "AUGGCGGAUGUUC"), type = "rna"), 1:13)
sq_rna_no_na_elem <- construct_sq(c("UAAACGGGCUA", "UACACGGGCGACU", "AGGCA", "AUGGCGGAUGUUC"), type = "rna")
sq_rna_no_na <- construct_sq(c("", "UACACGGGCGACU", "", "AUGGCGGAUGUUC"), type = "rna")

# PROTOTYPE PRESERVATION ----
test_that("remove_na() with only_elements = TRUE preserves all attributes of original vector", {
  expect_vector(remove_na(sq_ami, only_elements = TRUE),
                ptype = vec_ptype(sq_ami),
                size = vec_size(sq_ami))
  expect_vector(remove_na(sq_dna, only_elements = TRUE),
                ptype = vec_ptype(sq_dna),
                size = vec_size(sq_dna))
  expect_vector(remove_na(sq_rna, only_elements = TRUE),
                ptype = vec_ptype(sq_rna),
                size = vec_size(sq_rna))
})

test_that("remove_na() with only_elements = FALSE preserves all attributes of original vector", {
  expect_vector(remove_na(sq_ami, only_elements = FALSE),
                ptype = vec_ptype(sq_ami),
                size = vec_size(sq_ami))
  expect_vector(remove_na(sq_dna, only_elements = FALSE),
                ptype = vec_ptype(sq_dna),
                size = vec_size(sq_dna))
  expect_vector(remove_na(sq_rna, only_elements = FALSE),
                ptype = vec_ptype(sq_rna),
                size = vec_size(sq_rna))
})

# VALUE COMPUTATION ----
test_that("remove_na() with only_elements = TRUE returns correct value", {
  expect_equivalent(as.character(remove_na(sq_ami, only_elements = TRUE)),
                    as.character(sq_ami_no_na_elem))
  expect_equivalent(as.character(remove_na(sq_dna, only_elements = TRUE)),
                    as.character(sq_dna_no_na_elem))
  expect_equivalent(as.character(remove_na(sq_rna, only_elements = TRUE)),
                    as.character(sq_rna_no_na_elem))
})

test_that("remove_na() with only_elements = FALSE returns correct value", {
  expect_equivalent(as.character(remove_na(sq_ami, only_elements = FALSE)),
                    as.character(sq_ami_no_na))
  expect_equivalent(as.character(remove_na(sq_dna, only_elements = FALSE)),
                    as.character(sq_dna_no_na))
  expect_equivalent(as.character(remove_na(sq_rna, only_elements = FALSE)),
                    as.character(sq_rna_no_na))
})

# NO EFFECT ON SEQUENCES WITHOUT NA'S ----
test_that("remove_na() has no effect when original value has no NA's", {
  expect_identical(remove_na(sq_ami_no_na_elem, only_elements = TRUE),
                   sq_ami_no_na_elem)
  expect_identical(remove_na(sq_ami_no_na_elem, only_elements = FALSE),
                   sq_ami_no_na_elem)
})
