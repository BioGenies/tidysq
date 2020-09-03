sq_ami_cln <- clean(construct_sq(c("MIAANYTWIL","TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"),
                                 type = "ami"), only_elements = TRUE)

sq_ami_cln_short <- clean(construct_sq(c("MI","TI", "NY", "MA"), type = "ami"), only_elements = TRUE)
sq_ami_cln_short_2 <- clean(construct_sq(c("IAA","IAA", "YER", "AYI"), type = "ami"), only_elements = TRUE)

sq_ami_cln_cut <- clean(construct_sq(c("IAANYTWIL","IAALGNIIYRAIE", "YERTGHLI", "AYXXXIALN"),
                                     type = "ami"), only_elements = TRUE)

sq_ami_cln_cut_2 <- clean(construct_sq(c("ANYTWIL","ALGNIIYRAIE", "RTGHLI", "IALN"),
                                       type = "ami"), only_elements = TRUE)

sq_ami_cln_na <- structure(list(structure(as.raw(c(0xff, 0xff, 0xff, 0xff, 0xff)), original_length = 8),
                                structure(as.raw(c(0x88, 0xfc, 0xff, 0xff, 0xff)), original_length = 8),
                                structure(as.raw(c(0xff, 0xff, 0xff, 0xff, 0xff)), original_length = 8),
                                structure(as.raw(c(0xff, 0xff, 0xff, 0xff, 0xff)), original_length = 8)), 
                           class = c("clnsq", "amisq", "sq", "list"), 
                           alphabet = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-", "*"))

sq_ami_cln_na_2 <- structure(list(structure(as.raw(c(0xff, 0x7f)), original_length = 3),
                                  structure(as.raw(c(0x01, 0x11)), original_length = 3),
                                  structure(as.raw(c(0xff, 0x7f)), original_length = 3),
                                  structure(as.raw(c(0xff, 0x7f)), original_length = 3)), 
                             class = c("clnsq", "amisq", "sq", "list"), 
                             alphabet = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-", "*"))

#############

test_that("error on non-numeric 'indices' argument", {
  expect_error(bite(sq_ami_cln, "n"))
  expect_error(bite(sq_ami_cln, NA))
})

test_that("biting sq with positive indices - no reaching outside range", {
  expect_equal(bite(sq_ami_cln, 1:2),
               sq_ami_cln_short)
  expect_equal(bite(sq_ami_cln, 2:4),
               sq_ami_cln_short_2)
})

test_that("biting sq with positive indices - reaching outside range", {
  expect_warning(out <- bite(sq_ami_cln, 13:20), 
                 "some sequences are subsetted with index bigger than length - NA introduced")
  expect_equal(out, 
               sq_ami_cln_na)
  expect_warning(out <- bite(sq_ami_cln, 12:14), 
                 "some sequences are subsetted with index bigger than length - NA introduced")
  expect_equal(out, 
               sq_ami_cln_na_2)
})

test_that("biting sq with negative indices - no reaching outside range", {
  expect_equal(bite(sq_ami_cln, -1),
               sq_ami_cln_cut)
  expect_equal(bite(sq_ami_cln, -1:-3),
               sq_ami_cln_cut_2)
})

test_that("biting sq with negative indices - reaching outside range", {
  expect_equal(bite(sq_ami_cln, -15:-20),
               sq_ami_cln)
  expect_equal(bite(sq_ami_cln, -17:-18),
               sq_ami_cln)
})
