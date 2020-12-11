# SETUP ----
sq_dna <- sq(c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT"),
             alphabet = "dna_bsc")
sq_rna <- sq(c("", "KBS-UVW-AWWWG", "YGHHH-", "-CRASH", "MND-KUUBV-MY-"),
             alphabet = "rna_ext")
sq_ami <- sq(c("OUTLANDISH", "UNSTRUCTURIZED", "FEAR"),
             alphabet = "ami_ext")

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

# ERROR FOR UNSUPPORTED TARGET FORMAT ----
test_that("export_sq() throws an error whenever exporting to nonexistent format", {
  expect_error(export_sq(sq_dna, "MGMT::Kids"))
  expect_error(export_sq(sq_rna, "Darude::Sandstorm"))
  expect_error(export_sq(sq_ami, "Ferrari::F430"))
})
