# SETUP ----
name <- paste0("seq", 1:4)

str_dna_bsc <- c("CTAGTAG", "GGTAGATAG", "AAAAA")
str_ami_ext <- c("AMUILYBXX", "LCAMABAA", "F", "LVCGGA")
str_unt <- c("A4HH_AAX", "AHH1CPP")
str_atp <- c("mAmYmY", "nbAnsAmA")
long_sq <- sq(strrep("A", 1000), "dna_bsc")

alph_atp <- c("mA", "mY", "nbA", "nsA")

fasta_file_dna_bsc <- withr::local_tempfile()
writeLines(as.character(rbind(paste0(">", name[1:3]), str_dna_bsc)), 
           fasta_file_dna_bsc)

fasta_file_ami_ext <- withr::local_tempfile()
writeLines(as.character(rbind(paste0(">", name[1:4]), str_ami_ext)), 
           fasta_file_ami_ext)

fasta_file_unt <- withr::local_tempfile()
writeLines(as.character(rbind(paste0(">", name[1:2]), str_unt)), 
           fasta_file_unt)

fasta_file_atp <- withr::local_tempfile()
writeLines(as.character(rbind(paste0(">", name[1:2]), str_atp)), 
           fasta_file_atp)

fasta_file_blank_lines <- withr::local_tempfile()
writeLines(c("", ">sequence", "AGATA", "", "", ">sequence", "", "", "GAGAT"),
           fasta_file_blank_lines)

fasta_file_multiple_lines <- withr::local_tempfile()
writeLines(c(">sequence", "A", "C", "T", "G"),
           fasta_file_multiple_lines)

fasta_file_NA <- withr::local_tempfile()
writeLines(c(">sequence", "!!AC!!T!G!!A"),
           fasta_file_NA)

fasta_file_mixed_case <- withr::local_tempfile()
writeLines(c(">sequence", "aCTAgAGAAATGagATGAgAGGAT"),
           fasta_file_mixed_case)

# READING----
test_that("read_fasta() returns proper format of data", {
  fasta_dna <- read_fasta(fasta_file_dna_bsc, "dna_bsc")
  expect_s3_class(fasta_dna, "tbl_df", exact = FALSE)
  expect_s3_class(fasta_dna[["sq"]], "sq_dna_bsc", exact = FALSE)
  expect_true("character" %in% class(fasta_dna[["name"]]))
})

test_that("read_fasta() reads correct number of sequences", {
  expect_equal(nrow(read_fasta(fasta_file_dna_bsc, "dna_bsc")), 3)
  expect_equal(nrow(read_fasta(fasta_file_ami_ext, "ami_ext")), 4)
  expect_equal(nrow(read_fasta(fasta_file_unt, "unt")), 2)
  expect_equal(nrow(read_fasta(fasta_file_atp, alph_atp)), 2)
})


test_that("read_fasta() reads correctly sequences", {
  expect_equal(read_fasta(fasta_file_dna_bsc, "dna_bsc")[["sq"]],
               sq(str_dna_bsc, "dna_bsc"))
  expect_equal(read_fasta(fasta_file_ami_ext, "ami_ext")[["sq"]],
               sq(str_ami_ext, "ami_ext"))
  expect_equal(read_fasta(fasta_file_unt, "unt")[["sq"]],
               sq(str_unt, "unt"))
})

test_that("read_fasta() reads correctly name", {
  expect_equal(read_fasta(fasta_file_dna_bsc, "dna_bsc")[["name"]], name[1:3])
  expect_equal(read_fasta(fasta_file_ami_ext, "ami_ext")[["name"]], name[1:4])
  expect_equal(read_fasta(fasta_file_unt, "unt")[["name"]], name[1:2])
  expect_equal(read_fasta(fasta_file_atp, alph_atp)[["name"]], name[1:2])
})

test_that("read_fasta() skips blank lines", {
  expect_equal(read_fasta(fasta_file_blank_lines, "dna_bsc")[["sq"]],
               sq(c("AGATA", "GAGAT"), "dna_bsc"))
})

test_that("read_fasta() reads sequences with multiple lines", {
  expect_equal(read_fasta(fasta_file_multiple_lines, "dna_bsc")[["sq"]],
               sq("ACTG", "dna_bsc"))
})

test_that("read_fasta() detects type correctly", {
  expect_equal(read_fasta(fasta_file_dna_bsc, "dna_bsc"),
               read_fasta(fasta_file_dna_bsc))
  expect_equal(read_fasta(fasta_file_ami_ext, "ami_ext"),
               read_fasta(fasta_file_ami_ext))
})

test_that("read_fasta() reads NA values", {
  expect_equal(as.character(read_fasta(fasta_file_NA, "dna_bsc")[["sq"]]),
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
  expect_equal(read_fasta(fasta_file_atp, alph_atp)[["sq"]], 
               sq(str_atp, alph_atp))
})

# WRITING ----
test_that("write_fasta() creates a file at specified path", {
  withr::with_tempfile("fasta_out", {
    write_fasta(sq(str_dna_bsc, "dna_bsc"), name[1:3], fasta_out)
    expect_true(file.exists(fasta_out))
  })
})

test_that("write_fasta() saves sequences correctly", {
  withr::with_tempfile("fasta_out", {
    sq_dna <- sq(str_dna_bsc, "dna_bsc") 
    write_fasta(sq_dna, name[1:3], fasta_out)
    expect_equal(read_fasta(fasta_out, "dna_bsc")[["sq"]], sq_dna)
  })
})

test_that("write_fasta() saves names correctly", {
  withr::with_tempfile("fasta_out", {
    write_fasta(sq(str_ami_ext, "ami_ext"), name, fasta_out)
    expect_equal(read_fasta(fasta_out, "ami_ext")[["name"]], name)
  })
})

test_that("write_fasta() keeps line width", {
  withr::with_tempfile("fasta_out", {
    write_fasta(long_sq, name[1], fasta_out, width = 80)
    expect_true(all(nchar(readLines(fasta_out)) <= 80))
  })
})

# WRITE DATA.FRAME ----
test_that("data.frame columns are extracted and passed to write_fasta.sq()", {
  fasta_out_df <- withr::local_tempfile()
  write_fasta(read_fasta(fasta_file_dna_bsc, "dna_bsc"), fasta_out_df)
  
  fasta_out_sq <- withr::local_tempfile()
  write_fasta(sq(str_dna_bsc, "dna_bsc"), name[1:3], fasta_out_sq)
  
  expect_equal(
    read_fasta(fasta_out_df),
    read_fasta(fasta_out_sq)
  )
})

test_that("used data.frame columns can be specified", {
  fasta_out_df <- withr::local_tempfile()
  df_sq <- read_fasta(fasta_file_dna_bsc, "dna_bsc")
  df_sq[["name_upper"]] <- toupper(df_sq[["name"]])
  write_fasta(df_sq, fasta_out_df, .sq = "sq", .name = "name_upper")
  
  fasta_out_sq <- withr::local_tempfile()
  write_fasta(sq(str_dna_bsc, "dna_bsc"), toupper(name[1:3]), fasta_out_sq)
  
  expect_equal(
    read_fasta(fasta_out_df),
    read_fasta(fasta_out_sq)
  )
})

test_that("data.frame method properly passes 'width' parameter", {
  withr::with_tempfile("fasta_out", {
    write_fasta(
      data.frame(sq = long_sq, name = name[1]),
      fasta_out, width = 80
    )
    expect_true(all(nchar(readLines(fasta_out)) <= 80))
  })
})
