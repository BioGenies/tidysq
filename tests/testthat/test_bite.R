sq_ami_cln <- clean(construct_sq(c("MIAANYTWIL","TIAALGNIIYRAIE", 
                        "NYERTGHLI", "MAYXXXIALN"), type = "ami"), only_elements = TRUE)

sq_ami_cln_short <- clean(construct_sq(c("MI","TI", "NY", "MA"), type = "ami"), only_elements = TRUE)

sq_ami_cln_cut <- clean(construct_sq(c("IAANYTWIL","IAALGNIIYRAIE", 
                                       "YERTGHLI", "AYXXXIALN"), type = "ami"), only_elements = TRUE)

sq_ami_cln_na <- structure(list(as.raw(c(0xff, 0xff, 0xff, 0xff, 0xff)), 
                                as.raw(c(0x88, 0xfc, 0xff, 0xff, 0xff)), 
                                as.raw(c(0xff, 0xff, 0xff, 0xff, 0xff)), 
                                as.raw(c(0xff, 0xff, 0xff, 0xff, 0xff))), 
                           class = c("clnsq", "amisq", "sq"), 
                           alphabet = c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-", "*"))

#############

test_that("error on non-numeric 'indices' argument", 
          expect_error(bite(sq_ami_cln, "n")))

test_that("biting sq with positive indices - no reaching outside range", 
          expect_equal(bite(sq_ami_cln, 1:2),
                       sq_ami_cln_short))

test_that("biting sq with positive indices - reaching outside range", 
          expect_equal(bite(sq_ami_cln, 13:20), 
                       sq_ami_cln_na))

test_that("biting sq with negative indices - no reaching outside range", 
          expect_equal(bite(sq_ami_cln, -1),
                       sq_ami_cln_cut))

test_that("biting sq with negative indices - reaching outside range", 
          expect_equal(bite(sq_ami_cln, -15:-20),
                       sq_ami_cln))

