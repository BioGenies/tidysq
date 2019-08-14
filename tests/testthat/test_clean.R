sq_ami <- construct_sq(c("MIAANYTWIL","TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"), type = "ami")
sq_ami_cln <- construct_sq(c("MIAANYTWIL","TIAALGNIIYRAIE", "NYERTGHLI", ""), type = "ami")
sq_ami_cln_elements <- construct_sq(c("MIAANYTWIL","TIAALGNIIYRAIE", "NYERTGHLI", "MAYIALN"), type = "ami")
sq_nuc <- construct_sq(c("ATGCAGGA", "GACCGAACGAN", "TGACGAGCTTA", "ACTNNAGCN"), type = "nuc")
sq_nuc_cln <- construct_sq(c("ATGCAGGA", "", "TGACGAGCTTA", ""), type = "nuc")
sq_nuc_cln_elements <- construct_sq(c("ATGCAGGA", "GACCGAACGA", "TGACGAGCTTA", "ACTAGC"), type = "nuc")

test_that("cleaning doesn't affect clean ami sq", {
  expect_equal(clean(sq_ami_cln),
               sq_ami_cln)
})
          
test_that("cleaning doesn't affect clean nuc sq", {
  expect_equal(clean(sq_nuc_cln),
               sq_nuc_cln)
})
          
test_that("cleaning ami sq with parameter only_elements = FALSE", {
  expect_equal(clean(sq_ami, only_elements = FALSE),
               sq_ami_cln)
})
          
test_that("cleaning nuc sq with parameter only_elements = TRUE", {
  expect_equal(clean(sq_ami, only_elements = TRUE),
               sq_ami_cln_elements)
})
          
test_that("cleaning ami sq with parameter only_elements = FALSE", {
  expect_equal(clean(sq_nuc, only_elements = FALSE),
               sq_nuc_cln)
})
          
test_that("cleaning nuc sq with parameter only_elements = TRUE", {
  expect_equal(clean(sq_nuc, only_elements = TRUE),
               sq_nuc_cln_elements)
})
          
