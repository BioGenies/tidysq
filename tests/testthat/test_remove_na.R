sq_ami <- bite(construct_sq(c("MIAANYTWIL", "", "TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"), type = "ami"), 1:14)
sq_ami_removed_na_elements <- construct_sq(c("MIAANYTWIL", "", "TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"))
sq_ami_removed_na <- construct_sq(c("", "", "TIAALGNIIYRAIE", "", ""))
sq_dna <- bite(construct_sq(c("ATGCAGGA", "GACCGAACGAN", "TGACGAGCTTA", "ACTNNAGCN"), type = "dna"), 1:11)
sq_dna_removed_na_elements <- construct_sq(c("ATGCAGGA", "GACCGAACGAN", "TGACGAGCTTA", "ACTNNAGCN"), type = "dna")
sq_dna_removed_na <- construct_sq(c("", "", "TGACGAGCTTA", ""), type = "dna")
sq_rna <- bite(construct_sq(c("UAAACGGGCUA", "UACACGGGCGACU", "AGGCA", "AUGGCGGAUGUUC"), type = "rna"), 1:13)
sq_rna_removed_na_elements <- construct_sq(c("UAAACGGGCUA", "UACACGGGCGACU", "AGGCA", "AUGGCGGAUGUUC"), type = "rna")
sq_rna_removed_na <- construct_sq(c("", "UACACGGGCGACU", "", "AUGGCGGAUGUUC"), type = "rna")

test_that("removing NA with only_elements = TRUE", {
  expect_equal(remove_na(sq_ami, only_elements = TRUE),
               sq_ami_removed_na_elements)
  expect_equal(remove_na(sq_dna, only_elements = TRUE),
               sq_dna_removed_na_elements)
  expect_equal(remove_na(sq_rna, only_elements = TRUE),
               sq_rna_removed_na_elements)
})

test_that("removing NA with only_elements = FALSE", {
  expect_equal(clean(remove_na(sq_ami, only_elements = FALSE)),
               sq_ami_removed_na)
  expect_equal(clean(remove_na(sq_dna, only_elements = FALSE)),
               sq_dna_removed_na)
  expect_equal(clean(remove_na(sq_rna, only_elements = FALSE)),
               sq_rna_removed_na)
})
