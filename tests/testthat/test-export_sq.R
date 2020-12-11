# ERROR FOR NON-STANDARD SQ OBJECTS ----
test_that("export_sq() throws an error whenever passed object of class other that standard sq classes", {
  expect_error(export_sq(1:7, "Biostrings::DNAStringSet"))
  expect_error(export_sq(LETTERS, "Biostrings::DNAStringSet"))
  expect_error(export_sq(list(mean, sum, sd), "Biostrings::DNAStringSet"))
  expect_error(export_sq(sq(c(")R#)#!Vawr9fy", "*V)RUgBa^%#!b]")),
                         "Biostrings::DNAStringSet"))
  expect_error(export_sq(sq(c("accmsce", "auprcacc"), alphabet = c("auprc", "acc", "msce")),
                         "Biostrings::DNAStringSet"))
})
