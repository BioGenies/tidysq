file_ami <- system.file(package = "tidysq", "/sample_fasta/sample_ami.fasta")
file_ami_nonst <- system.file(package = "tidysq", "/sample_fasta/sample_ami_nonst.fasta")
file_cln_ami <- system.file(package = "tidysq", "/sample_fasta/sample_cln_ami.fasta")
file_cln_dna <- system.file(package = "tidysq", "/sample_fasta/sample_cln_dna.fasta")
file_dna <- system.file(package = "tidysq", "/sample_fasta/sample_dna.fasta")

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
          
test_that("reading fasta file with clean DNA sequences with type and is_clean unspecified", {
  expect_silent(read_fasta(file_cln_dna))
  expect_equal(dim(read_fasta(file_cln_dna)),
               c(100, 2))
  expect_equal(read_fasta(file_cln_dna)[[2]][1],
               construct_sq("CTTTATGATTTCTTGTAGTACTCCTCTTGA", type = "dna"))
})
          
test_that("reading fasta file with ambiguous DNA sequences with type and is_clean unspecified", {
  expect_silent(read_fasta(file_dna))
  expect_equal(dim(read_fasta(file_dna)),
               c(100, 2))
  expect_equal(read_fasta(file_dna)[[2]][1],
               construct_sq("MBSDCKKHSGGMNGTYVKWKWGCVWTYM", type = "dna"))
})
          
test_that("reading fasta file with non-standard element", {
  expect_silent(read_fasta(file_ami_nonst))
  expect_equal(dim(read_fasta(file_ami_nonst)),
               c(100, 2))
  skip("cannot reconstruct original structure")
  expect_equal(read_fasta(file_ami_nonst)[[2]][1],
               construct_sq("NFGmAERESWIQIMSFIWHRVNNLAYQPQH", type = "unt"))
})
