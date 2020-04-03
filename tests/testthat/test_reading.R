file_ami <- system.file(package = "tidysq", "/sample_fasta/sample_ami.fasta")
file_ami_nonst <- system.file(package = "tidysq", "/sample_fasta/sample_ami_nonst.fasta")
file_cln_ami <- system.file(package = "tidysq", "/sample_fasta/sample_cln_ami.fasta")
file_cln_nuc <- system.file(package = "tidysq", "/sample_fasta/sample_cln_nuc.fasta")
file_nuc <- system.file(package = "tidysq", "/sample_fasta/sample_nuc.fasta")

test_that("reading fasta file with clean aminoacid sequences with type and is_clean unspecified", {
  expect_silent(read_fasta(file_cln_ami))
  expect_equal(dim(read_fasta(file_cln_ami)),
               c(100, 2))
  expect_equal(read_fasta(file_cln_ami)[[2]][1],
               construct_sq("FKFNDTEMQAHFEFHFKWTSFCCDTDGWGTN", type = "ami"))
})

test_that("reading fasta file with ambiguous aminoacid sequences with type and is_clean unspecified", {
  expect_silent(read_fasta(file_ami))
  expect_equal(dim(read_fasta(file_ami)),
               c(100, 2))
  expect_equal(read_fasta(file_ami)[[2]][1],
               construct_sq("KWOVPFWSAZFBTPSQBIKFBQDQXAFNY", type = "ami"))
})
          
test_that("reading fasta file with clean nucleotides sequences with type and is_clean unspecified", {
  expect_silent(read_fasta(file_cln_nuc))
  expect_equal(dim(read_fasta(file_cln_nuc)),
               c(100, 2))
  expect_equal(read_fasta(file_cln_nuc)[[2]][1],
               construct_sq("CTTTATGATTTCTTGTAGTACTCCUCTTGA", type = "nuc"))
})
          
test_that("reading fasta file with ambiguous nucleotides sequences with type and is_clean unspecified", {
  expect_silent(read_fasta(file_nuc))
  expect_equal(dim(read_fasta(file_nuc)),
               c(100, 2))
  expect_equal(read_fasta(file_nuc)[[2]][1],
               construct_sq("MBSDCKKHSGGMNGTYVKWKWGCVWUYM", type = "nuc"))
})
          
test_that("reading fasta file with non-standard element", {
  expect_silent(read_fasta(file_ami_nonst))
  expect_equal(dim(read_fasta(file_ami_nonst)),
               c(100, 2))
  skip("cannot reconstruct original structure")
  expect_equal(read_fasta(file_ami_nonst)[[2]][1],
               construct_sq("NFGmAERESWIQIMSFIWHRVNNLAYQPQH", type = "unt"))
})
