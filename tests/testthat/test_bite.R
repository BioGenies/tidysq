sq_ami_cln <- clean(construct_sq(c("MIAANYTWIL","TIAALGNIIYRAIE", 
                        "NYERTGHLI", "MAYXXXIALN"), type = "ami"), only_elements = TRUE)

sq_ami_cln_short <- clean(construct_sq(c("MI","TI", "NY", "MA"), type = "ami"), only_elements = TRUE)

sq_ami_cln_cut <- clean(construct_sq(c("IAANYTWIL","IAALGNIIYRAIE", 
                                       "YERTGHLI", "AYXXXIALN"), type = "ami"), only_elements = TRUE)

##

test_that("error on non-numeric 'indices' argument", 
          expect_error(bite(sq_ami_cln, "n")))

test_that("biting sq with positive indices - no reaching outside range", 
          expect_equal(bite(sq_ami_cln, 1:2),
                       sq_ami_cln_short))

test_that("biting sq with positive indices - reaching outside range", {})

test_that("biting sq with negative indices - no reaching outside range", 
          expect_equal(bite(sq_ami_cln, -1),
                       sq_ami_cln_cut))

test_that("biting sq with negative indices - reaching outside range", {})

