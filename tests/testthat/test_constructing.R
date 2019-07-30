test_that("single short ambiguous nucleotides sequences constructing", {
  expect_identical(construct_sq("ATGN"),
                   structure(list(as.raw(c(0x81, 0x0c, 0x08))), 
                             alphabet = c("A", "C", "G", "T", "U", "W", 
                                          "S", "M", "K", "R", "Y", "B", 
                                          "D", "H", "V", "N", "-"), 
                             class = c("nucsq", "sq")))
})
test_that("multiple short ambiguous nucleotides sequences constructing", {})
test_that("single medium-length ambiguous nucleotides sequences constructing", {})
test_that("multiple medium-length ambiguous nucleotides sequences constructing", {})
test_that("single long ambiguous nucleotides sequences constructing", {})
test_that("multiple long ambiguous nucleotides sequences constructing", {})
test_that("single short non-ambiguous nucleotides sequences constructing", {
  expect_identical(construct_sq("ATG"),
                   structure(list(as.raw(c(0xe1, 0x00))), 
                             alphabet = c("A", "C", "G", "T", "U", "-"), 
                             class = c("clnsq", "nucsq", "sq")))
})
test_that("multiple short non-ambiguous nucleotides sequences constructing", {})
test_that("single medium-length non-ambiguous nucleotides sequences constructing", {})
test_that("multiple medium-length non-ambiguous nucleotides sequences constructing", {
  expect_identical(construct_sq(c("CTGATGGATGATGCTAGC", 
                                  "TCTGATGATGAAGCTAGTAGGAA", 
                                  "GGGTAGGATAGGCTAGATAG")),
                   structure(list(as.raw(c(0xe2, 0xc2, 0x2d, 0x5c, 0x38, 0x31, 0x13)), 
                                  as.raw(c(0x14, 0x17, 0x2e, 0x5c, 0x32, 0x31, 0x63, 0xb6, 0x04)), 
                                  as.raw(c(0xdb, 0x98, 0x2d, 0xcc, 0x26, 0x66, 0x61, 0x06))), 
                             alphabet = c("A", "C", "G", "T", "U", "-"), 
                             class = c("clnsq", "nucsq", "sq")))
})
test_that("single long non-ambiguous nucleotides sequences constructing", {})
test_that("multiple long non-ambiguous nucleotides sequences constructing", {})
test_that("single short ambiguous aminoacids sequences constructing", {
  expect_identical(construct_sq("ATGN"),
                   structure(list(as.raw(c(0x81, 0x0c, 0x08))), 
                             alphabet = c("A", "C", "G", "T", "U", "W", 
                                          "S", "M", "K", "R", "Y", "B", 
                                          "D", "H", "V", "N", "-"), 
                             class = c("nucsq", "sq")))
})
test_that("multiple short ambiguous aminoacids sequences constructing", {})
test_that("single medium-length ambiguous aminoacids sequences constructing", {})
test_that("multiple medium-length ambiguous aminoacids sequences constructing", {})
test_that("single long ambiguous aminoacids sequences constructing", {})
test_that("multiple long ambiguous aminoacids sequences constructing", {})
test_that("single short non-ambiguous aminoacids sequences constructing", {
  expect_identical(construct_sq("PQKH"),
                   structure(list(as.raw(c(0xcd, 0xa5, 0x03))), 
                             alphabet = c("A", "C", "D", "E", "F", "G", 
                                          "H", "I", "K", "L", "M", "N", 
                                          "P", "Q", "R", "S", "T", "V", 
                                          "W", "Y", "-"), 
                             class = c("clnsq", "amisq", "sq")))
})
test_that("multiple short non-ambiguous aminoacids sequences constructing", {})
test_that("single medium-length non-ambiguous aminoacids sequences constructing", {})
test_that("multiple medium-length non-ambiguous aminoacids sequences constructing", {})
test_that("single long non-ambiguous aminoacids sequences constructing", {})
test_that("multiple long non-ambiguous aminoacids sequences constructing", {})

test_that("ambiguous nucleotides with specified type sequences constructing", {
  expect_equal(construct_sq("GGNWGATCGANN", type = "nuc"),
               structure(list(as.raw(c(0x63, 0x40, 0x33, 0x02, 0x11, 0x23, 0x40, 0x08))), 
                         alphabet = c("A", "C", "G", "T", "U", "W", 
                                      "S", "M", "K", "R", "Y", "B", 
                                      "D", "H", "V", "N", "-"), 
                         class = c("nucsq", "sq")))
})
test_that("non-ambiguous nucleotides with specified type sequences constructing", {})
test_that("ambiguous aminoacids with specified type sequences constructing", {})
test_that("non-ambiguous aminoacids with specified type sequences constructing", {})
test_that("ambiguous nucleotides with specified is_clean sequences constructing", {})
test_that("non-ambiguous nucleotides with specified is_clean sequences constructing", {})
test_that("ambiguous aminoacids with specified is_clean sequences constructing", {})
test_that("non-ambiguous aminoacids with specified is_clean sequences constructing", {})
test_that("ambiguous nucleotides with specified both type and is_clean sequences constructing", {})
test_that("non-ambiguous nucleotides with specified both type and is_clean sequences constructing", {})
test_that("ambiguous aminoacids with specified both type and is_clean sequences constructing", {})
test_that("non-ambiguous aminoacids with specified both type and is_clean sequences constructing", {})