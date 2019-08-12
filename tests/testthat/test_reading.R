file_ami <- system.file(package = "tidysq", "/sample_fasta/sample_ami.fasta")
file_ami_nonst <- system.file(package = "tidysq", "/sample_fasta/sample_ami_nonst.fasta")
file_cln_ami <- system.file(package = "tidysq", "/sample_fasta/sample_cln_ami.fasta")
file_cln_nuc <- system.file(package = "tidysq", "/sample_fasta/sample_cln_nuc.fasta")
file_nuc <- system.file(package = "tidysq", "/sample_fasta/sample_nuc.fasta")

test_that("reading fasta file with clean aminoacid sequences with type and is_clean unspecified", 
          expect_silent(read_fasta(file_cln_ami)))
test_that("reading fasta file with ambiguous aminoacid sequences with type and is_clean unspecified", 
          expect_silent(read_fasta(file_ami)))
test_that("reading fasta file with clean nucleotides sequences with type and is_clean unspecified", 
          expect_silent(read_fasta(file_cln_nuc)))
test_that("reading fasta file with ambiguous nucleotides sequences with type and is_clean unspecified", 
          expect_silent(read_fasta(file_nuc)))

test_that("reading fasta file with non-standard element", 
          expect_silent(read_fasta(file_ami_nonst)))
