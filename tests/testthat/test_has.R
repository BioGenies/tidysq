sq_ami <- construct_sq(c("AGNTYIKFGGAYTIB", "MATEGILIAADGYTWIL", "MIPADHICAANGIENAGIK"), type = 'ami')
sq_nuc <- construct_sq(c("CTGAATGCAGTACCGTAAT", "ATGCCGTAAATGCCAT", "CAGACCANNNATAG"), type = 'nuc')

test_that("%has% detects correctly motif that is single unambiguous amino acid in sequences", {
  expect_equal(sum(sq_ami %has% "A"),
               3)
  expect_equal(sum(sq_ami %has% "C"),
               1)
})

test_that("%has% detects correctly motif that is single unambiguous nucleotide in sequences", {
  expect_equal(sum(sq_nuc %has% "A"),
               3)
  expect_equal(sum(sq_nuc %has% "C"),
               3)
})

test_that("%has% detects correctly motif that is single ambiguous amino acid in sequences", {
  expect_equal(sum(sq_ami %has% "B"), 
               3)
  expect_equal(sum(sq_ami %has% "J"), 
               3)
})

test_that("%has% detects correctly motif that is single ambiguous nucleotide in sequences", {
  expect_equal(sum(sq_nuc %has% "W"),
               3)
  expect_equal(sum(sq_nuc %has% "K"),
               3)
})

test_that("%has% detects correctly leading letters of amino acid sequences using '^'", {
  expect_equal(sum(sq_ami %has% "^A"),
               1)
  expect_equal(sum(sq_ami %has% "^M"),
               2)
})

test_that("%has% detects correctly leading letters of nucleotide sequences using '^'", {
  expect_equal(sum(sq_nuc %has% "^C"), 
               2)
  expect_equal(sum(sq_nuc %has% "^A"), 
               1)
})

test_that("%has% detects correctly letters at the end of amino acid sequences using '$'", {
  expect_equal(sum(sq_ami %has% "K$"), 
               1)
  expect_equal(sum(sq_ami %has% "L$"), 
               1)
  expect_equal(sum(sq_ami %has% "B$"), 
               1)
})
          
test_that("%has% detects correctly letters at the end of nucleotide sequences using '$'", {
  expect_equal(sum(sq_nuc %has% "T$"),
               2)
  expect_equal(sum(sq_nuc %has% "G$"),
               1)
})
          
test_that("%has% detects correctly multiple-letter motifs in amino acid sequences", {
  expect_equal(sum(sq_ami %has% "TY"),
               1)
  expect_equal(sum(sq_ami %has% "GAY"),
               1)
  expect_equal(sum(sq_ami %has% "MATE"),
               1)
})

test_that("%has% detects correctly multiple-letter motifs in nucleotide sequences", {
  expect_equal(sum(sq_nuc %has% "CC"),
               3)
  expect_equal(sum(sq_nuc %has% "TAA"),
               2)
  expect_equal(sum(sq_nuc %has% "TAAT"),
               1)
})

test_that("%has% detects correctly multiple motifs in amino acid sequences", {
  expect_equal(sum(sq_ami %has% c("TY", "GA")),
               1)
  expect_equal(sum(sq_ami %has% c("MATE", "TY")),
               0)
  expect_equal(sum(sq_ami %has% c("GI", "MATE")),
               1)
})
          
test_that("%has% detects correctly multiple motifs in nucleotide sequences", {
  expect_equal(sum(sq_nuc %has% c("CC", "AT")),
               3)
  expect_equal(sum(sq_nuc %has% c("CC", "TAAT")),
               1)
  expect_equal(sum(sq_nuc %has% c("CC", "TAA")),
               2)
})
        