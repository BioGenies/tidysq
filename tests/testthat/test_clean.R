sq_ami <- construct_sq(c("MIAANYTWIL","TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"), type = "ami")
sq_ami_cln <- construct_sq(c("MIAANYTWIL","TIAALGNIIYRAIE", "NYERTGHLI", ""), type = "ami")
sq_ami_cln_elements <- construct_sq(c("MIAANYTWIL","TIAALGNIIYRAIE", "NYERTGHLI", "MAYIALN"), type = "ami")
sq_dna <- construct_sq(c("ATGCAGGA", "GACCGAACGAN", "TGACGAGCTTA", "ACTNNAGCN"), type = "dna")
sq_dna_cln <- construct_sq(c("ATGCAGGA", "", "TGACGAGCTTA", ""), type = "dna")
sq_dna_cln_elements <- construct_sq(c("ATGCAGGA", "GACCGAACGA", "TGACGAGCTTA", "ACTAGC"), type = "dna")
sq_rna <- construct_sq(c("UGGCGNNGBV", "ACGGUUUCGUU", "UBBGDGAACG", "GGCUCGACAGACUGC"), type = "rna")
sq_rna_cln <- construct_sq(c("", "ACGGUUUCGUU", "", "GGCUCGACAGACUGC"), type = "rna")
sq_rna_cln_elements <- construct_sq(c("UGGCGG", "ACGGUUUCGUU", "UGGAACG", "GGCUCGACAGACUGC"), type = "rna")

test_that("cleaning doesn't affect clean ami sq", {
  expect_equal(clean(sq_ami_cln),
               sq_ami_cln)
})
test_that("cleaning doesn't affect clean dna sq", {
  expect_equal(clean(sq_dna_cln),
               sq_dna_cln)
})
test_that("cleaning doesn't affect clean rna sq", {
  expect_equal(clean(sq_rna_cln),
               sq_rna_cln)
})

test_that("cleaning ami sq with parameter only_elements = FALSE", {
  expect_equal(clean(sq_ami, only_elements = FALSE),
               sq_ami_cln)
})
test_that("cleaning dna sq with parameter only_elements = FALSE", {
  expect_equal(clean(sq_dna, only_elements = FALSE),
               sq_dna_cln)
})
test_that("cleaning rna sq with parameter only_elements = FALSE", {
  expect_equal(clean(sq_rna, only_elements = FALSE),
               sq_rna_cln)
})

test_that("cleaning ami sq with parameter only_elements = TRUE", {
  expect_equal(clean(sq_ami, only_elements = TRUE),
               sq_ami_cln_elements)
})
test_that("cleaning dna sq with parameter only_elements = TRUE", {
  expect_equal(clean(sq_dna, only_elements = TRUE),
               sq_dna_cln_elements)
})
test_that("cleaning rna sq with parameter only_elements = TRUE", {
  expect_equal(clean(sq_rna, only_elements = TRUE),
               sq_rna_cln_elements)
})
