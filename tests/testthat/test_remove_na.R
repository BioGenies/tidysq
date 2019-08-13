sq_ami <- bite(construct_sq(c("MIAANYTWIL", "", "TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN")), 1:14)
sq_ami_removed_na_elements <- construct_sq(c("MIAANYTWIL", "", "TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"))
sq_ami_removed_na <- construct_sq(c("", "", "TIAALGNIIYRAIE", "", ""))

test_that("removing NA with only_elements = TRUE", 
          expect_equal(remove_na(sq_ami, only_elements = TRUE),
                       sq_ami_removed_na_elements))
test_that("removing NA with only_elements = FALSE", 
          expect_equal(clean(remove_na(sq_ami, only_elements = FALSE)),
                       sq_ami_removed_na))