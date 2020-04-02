library(tibble)

ami_sq <- construct_sq(c("ALAALGALAGLL", "GACNNALA", "AALGLA"))
name <- c("p1", "p2", "p3")

result <- tibble(name = c("p1", "p1", "p1", "p1", "p2", "p2", "p3", "p3"),
                 peptides = construct_sq(c("AL", "AAL", "GAL", "AGLL", "GACNNAL", "A", "AAL", "GLA")))

test_that("digesting sq with simple pattern", {
  expect_equal(digest(ami_sq, name, "[AL]."), result)
})