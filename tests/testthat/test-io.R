# SETUP ----
name <- paste0("seq", 1:4)

str_dna_bsc <- c("CTAGTAG", "GGTAGATAG", "AAAAA")
str_ami_ext <- c("AMUILYBXX", "LCAMABAA", "F", "LVCGGA")
str_unt <- c("A4HH_AAX", "AHH1CPP")
str_atp <- c("mAmYmY", "nbAnsAmA")
long_sq <- sq(strrep("A", 1000), "dna_bsc")

alph_atp <- c("mA", "mY", "nbA", "nsA")

fasta_file_dna_bsc <- tempfile()
writeLines(as.character(rbind(paste0(">", name[1:3]), str_dna_bsc)), 
           fasta_file_dna_bsc)

fasta_file_ami_ext <- tempfile()
writeLines(as.character(rbind(paste0(">", name[1:4]), str_ami_ext)), 
           fasta_file_ami_ext)

fasta_file_unt <- tempfile()
writeLines(as.character(rbind(paste0(">", name[1:2]), str_unt)), 
           fasta_file_unt)

fasta_file_atp <- tempfile()
writeLines(as.character(rbind(paste0(">", name[1:2]), str_atp)), 
           fasta_file_atp)


fasta_file_blank_lines <- tempfile()
writeLines(c("", ">sequence", "AGATA", "", "", ">sequence", "", "", "GAGAT"),
           fasta_file_blank_lines)

fasta_file_multiple_lines <- tempfile()
writeLines(c(">sequence", "A", "C", "T", "G"),
           fasta_file_multiple_lines)

fasta_file_NA <- tempfile()
writeLines(c(">sequence", "!!AC!!T!G!!A"),
           fasta_file_NA)

fasta_file_mixed_case <- tempfile()
writeLines(c(">sequence", "aCTAgAGAAATGagATGAgAGGAT"),
           fasta_file_mixed_case)

fasta_file_out_1 <- paste0(tempdir(), "/tidysq_fasta_out_1")
fasta_file_out_2 <- paste0(tempdir(), "/tidysq_fasta_out_2")
fasta_file_out_3 <- paste0(tempdir(), "/tidysq_fasta_out_3")
fasta_file_out_4 <- paste0(tempdir(), "/tidysq_fasta_out_4")

# READING----
test_that("read_fasta() returns proper format of data", {
  fasta_dna <- read_fasta(fasta_file_dna_bsc, "dna_bsc")
  expect_s3_class(fasta_dna, "tbl_df", exact = FALSE)
  expect_s3_class(fasta_dna$sq, "sq_dna_bsc", exact = FALSE)
  expect_true("character" %in% class(fasta_dna$name))
})

test_that("read_fasta() reads correct number of sequences", {
  expect_equal(nrow(read_fasta(fasta_file_dna_bsc, "dna_bsc")), 3)
  expect_equal(nrow(read_fasta(fasta_file_ami_ext, "ami_ext")), 4)
  expect_equal(nrow(read_fasta(fasta_file_unt, "unt")), 2)
  expect_equal(nrow(read_fasta(fasta_file_atp, alph_atp)), 2)
})


test_that("read_fasta() reads correctly sequences", {
  expect_equal(read_fasta(fasta_file_dna_bsc, "dna_bsc")$sq, 
               sq(str_dna_bsc, "dna_bsc"))
  expect_equal(read_fasta(fasta_file_ami_ext, "ami_ext")$sq, 
               sq(str_ami_ext, "ami_ext"))
  expect_equal(read_fasta(fasta_file_unt, "unt")$sq, 
               sq(str_unt, "unt"))
})

test_that("read_fasta() reads correctly name", {
  expect_equal(read_fasta(fasta_file_dna_bsc, "dna_bsc")$name, name[1:3])
  expect_equal(read_fasta(fasta_file_ami_ext, "ami_ext")$name, name[1:4])
  expect_equal(read_fasta(fasta_file_unt, "unt")$name, name[1:2])
  expect_equal(read_fasta(fasta_file_atp, alph_atp)$name, name[1:2])
})

test_that("read_fasta() skips blank lines", {
  expect_equal(read_fasta(fasta_file_blank_lines, "dna_bsc")$sq,
               sq(c("AGATA", "GAGAT"), "dna_bsc"))
})

test_that("read_fasta() read sequences with multiple lines", {
  expect_equal(read_fasta(fasta_file_multiple_lines, "dna_bsc")$sq,
               sq(c("ACTG"), "dna_bsc"))
})

test_that("read_fasta() detects correctly type", {
  expect_equal(read_fasta(fasta_file_dna_bsc, "dna_bsc"),
               read_fasta(fasta_file_dna_bsc))
  expect_equal(read_fasta(fasta_file_ami_ext, "ami_ext"),
               read_fasta(fasta_file_ami_ext))
})

test_that("read_fasta() reads NA values", {
  expect_equal(as.character(read_fasta(fasta_file_NA, "dna_bsc")$sq),
               "!!AC!!T!G!!A")
})

test_that("read_fasta() throws warning when ignore_case = FALSE and safe_mode = TRUE", {
  expect_warning(read_fasta(fasta_file_mixed_case, "dna_bsc", safe_mode = TRUE), 
                 "Detected letters that do not match specified type!")
})

test_that("read_fasta() throws warning when strange characters detected and safe_mode = TRUE", {
  expect_warning(read_fasta(fasta_file_unt, "dna_bsc", safe_mode = TRUE), 
                 "Detected letters that do not match specified type!")
})

test_that("read_fasta() reads multichar letters correctly", {
  expect_equal(read_fasta(fasta_file_atp, alph_atp)$sq, 
               sq(str_atp, alph_atp))
})

# WRITING ----

test_that("write_fasta() creates a file at specified path", {
  write_fasta(sq(str_dna_bsc, "dna_bsc"), 
              name[1:3], fasta_file_out_1)
  expect_true(file.exists(fasta_file_out_1))
})

test_that("write_fasta() saves sequences correctly", {
  sq_dna <- sq(str_dna_bsc, "dna_bsc") 
  write_fasta(sq_dna, 
              name[1:3], fasta_file_out_2)
  expect_equal(read_fasta(fasta_file_out_2, "dna_bsc")$sq, sq_dna)
})

test_that("write_fasta() saves names correctly", {
  write_fasta(sq(str_ami_ext, "ami_ext"), name, fasta_file_out_3)
  expect_equal(read_fasta(fasta_file_out_3, "ami_ext")$name, name)
})

test_that("write_fasta() keeps line width", {
  write_fasta(long_sq, name[1], fasta_file_out_4, width = 80)
  expect_true(all(nchar(readLines(fasta_file_out_4)) <= 80))
})

# CLEANUP ----

file.remove(fasta_file_dna_bsc,
            fasta_file_ami_ext,
            fasta_file_unt,
            fasta_file_atp,
            fasta_file_blank_lines,
            fasta_file_multiple_lines,
            fasta_file_NA,
            fasta_file_mixed_case,
            fasta_file_out_1,
            fasta_file_out_2,
            fasta_file_out_3,
            fasta_file_out_4)
