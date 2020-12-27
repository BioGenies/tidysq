# SETUP ----
sq_dna_bsc <- sq(c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT"),
                 alphabet = "dna_bsc")
sq_dna_ext <- sq(c("NARTYVTCY", "", "ATKCYGDD", "", "DNAKYTD"),
                 alphabet = "dna_ext")
sq_rna_ext <- sq(c("", "KBS-UVW-AWWWG", "YGHHH-", "-CRASH", "MND-KUUBV-MY-"),
                 alphabet = "rna_ext")
sq_ami_bsc <- sq(c("ACEH", "PASAI", "MALACCA", "SIAK"),
                 alphabet = "ami_bsc")
sq_unt <- sq(c("VIP01", "VIP02", "VIP04", "MISSING_ONE"),
             alphabet = "unt")
sq_atp <- sq(c("mAmYmY", "nbAnsAmA", ""),
             alphabet = c("mA", "mY", "nbA", "nsA"))

# Names below are just chosen for extravagance, not political reasons
names_5 <- c("Monza", "Imola", "Mugello", "Pescara", "Modena")
names_4 <- c("gara", "zara", "zarete", "dira")
names_3 <- c("naiz", "haiz", "da")

N_interpretations <- c("A", "C", "G", "T", "K", "R", "Y", "W", "S", "M", "B", "D", "H", "V", "N")

# ERROR FOR NON-SQ OBJECTS ----
test_that("find_motifs() throws an error whenever passed object of class other that sq", {
  expect_error(find_motifs(1:5, names_5, "ACC"))
  expect_error(find_motifs(LETTERS, letters, "H"))
  expect_error(find_motifs(list(mean, sum, sd), names_3, "mean"))
})

# CORRECT PROTOTYPE OF RETURNED VALUE ----
test_that("find_motifs() returns a tibble with columns specified in docs", {
  expect_vector(find_motifs(sq_dna_bsc, names_5, "TAG"),
                ptype = tibble::tibble(
                  names = character(),
                  found = vec_ptype(sq_dna_bsc),
                  sought = character(),
                  start = integer(),
                  end = integer()
                ))
  expect_vector(find_motifs(sq_rna_ext, names_5, c("YNG", "CHK")),
                ptype = tibble::tibble(
                  names = character(),
                  found = vec_ptype(sq_rna_ext),
                  sought = character(),
                  start = integer(),
                  end = integer()
                ))
  expect_vector(find_motifs(sq_unt, names_4, "^VIP"),
                ptype = tibble::tibble(
                  names = character(),
                  found = vec_ptype(sq_unt),
                  sought = character(),
                  start = integer(),
                  end = integer()
                ))
  expect_vector(find_motifs(sq_atp, names_3, "mYmY"),
                ptype = tibble::tibble(
                  names = character(),
                  found = vec_ptype(sq_atp),
                  sought = character(),
                  start = integer(),
                  end = integer()
                ))
})

# ARGUMENT PREREQUISITES ----
test_that("x must be an sq object", {
  expect_error(find_motifs(list(5, LETTERS, mean), names_3, "TAG"),
               "method 'find_motifs' isn't implemented for this type of object")
})

test_that("name argument must contain unique elements", {
  # There would be a test that the same code works with unique names, but
  # it would involve creating separate expectation and we don't have to do that,
  # as said code is already tested in prototype tests
  expect_error(find_motifs(sq_dna_bsc,
                           c("Monza", "Imola", "Mugello", "Monza", "Mugello"),
                           "TAG"))
})

test_that("^ cannot appear in any position other than first", {
  expect_error(find_motifs(sq_dna_bsc, names_5, "T^A"))
  expect_error(find_motifs(sq_unt, names_4, c("^VIP", "ONE^")))
})

test_that("$ cannot appear in any position other than last", {
  expect_error(find_motifs(sq_dna_bsc, names_5, "T$A"))
  expect_error(find_motifs(sq_unt, names_4, c("^VIP", "$ONE")))
})

# NAMES COLUMN ----
test_that("'names' column contains only elements from vector passed to find_motifs()", {
  expect_subset(find_motifs(sq_dna_bsc, names_5, "TAG")[["names"]],
                names_5)
  expect_subset(find_motifs(sq_rna_ext, names_5, c("YNG", "CHK"))[["names"]],
                names_5)
  expect_subset(find_motifs(sq_unt, names_4, "^VIP")[["names"]],
                names_4)
  expect_subset(find_motifs(sq_atp, names_3, "mYmY")[["names"]],
                names_3)
})

# SOUGHT COLUMN ----
test_that("'sought' column contains a subset of searched motifs", {
  expect_subset(find_motifs(sq_dna_bsc, names_5, "TAG")[["sought"]],
                "TAG")
  expect_subset(find_motifs(sq_dna_ext, names_5, "NND$")[["sought"]],
                "NND$")
  expect_subset(find_motifs(sq_ami_bsc, names_4, c("AC", "AI", "AH"))[["sought"]],
                c("AC", "AI", "AH"))
  expect_subset(find_motifs(sq_rna_ext, names_5, c("YNG", "CHK"))[["sought"]],
                c("YNG", "CHK"))
  expect_subset(find_motifs(sq_unt, names_4, "^VIP")[["sought"]],
                "^VIP")
  expect_subset(find_motifs(sq_atp, names_3, "mYmY")[["sought"]],
                "mYmY")
})

# FOUND COLUMN ----
test_that("'found' column contains sequences that are subset of searched simple motifs", {
  expect_subset(as.character(find_motifs(sq_dna_bsc, names_5, "TAG")[["found"]]),
                "TAG")
  expect_subset(as.character(find_motifs(sq_rna_ext, names_5, "AU")[["found"]]),
                "AU")
  expect_subset(as.character(find_motifs(sq_ami_bsc, names_4, c("AC", "AI", "AH"))[["found"]]),
                c("AC", "AI", "AH"))
  expect_subset(as.character(find_motifs(sq_unt, names_4, c("OUGH", "VIP"))[["found"]]),
                c("OUGH", "VIP"))
  expect_subset(as.character(find_motifs(sq_atp, names_3, c("mAmY", "nsA"))[["found"]]),
                c("mAmY", "nsA"))
})

test_that("special characters are not included in 'found' column", {
  expect_subset(as.character(find_motifs(sq_dna_bsc, names_5, "G$")[["found"]]),
                "G")
  expect_subset(as.character(find_motifs(sq_unt, names_4, c("^VIP", "ONE$", "_"))[["found"]]),
                c("VIP", "ONE", "_"))
  expect_subset(as.character(find_motifs(sq_atp, names_3, "mYmY$")[["found"]]),
                "mYmY")
})

test_that("'found' column can contain any interpretation of an ambiguous motif", {
  expect_subset(as.character(find_motifs(sq_dna_bsc, names_5, "ANT")[["found"]]),
                paste0("A", N_interpretations, "T"))
  expect_subset(as.character(find_motifs(sq_ami_bsc, names_4, c("XAI", "CX"))[["found"]]),
                c(paste0(LETTERS, "AI"), paste0("C", LETTERS)))
})

test_that("'found' column can handle both special characters and ambiguous letters at once", {
  expect_subset(as.character(find_motifs(sq_ami_bsc, names_4, "^XA")[["found"]]),
                paste0(LETTERS, "A"))
})

# INDEX COLUMNS ----
test_that("'start' and 'end' columns have values between 1 and length(sequence)", {
  sqibble_1 <- find_motifs(sq_dna_bsc, names_5, "TAG")
  sqibble_1[["found_length"]] <- get_sq_lengths(sqibble_1[["found"]])
  purrr::pwalk(sqibble_1, function(names, sought, found, start, end, found_length) {
    sequence_length <- get_sq_lengths(sq_dna_bsc)[which(names == names_5)]
    expect_gte(start, 1)
    expect_gte(end - found_length + 1, 1)
    expect_lte(start + found_length - sequence_length, 1)
    expect_lte(end - sequence_length + 1, 1)
  })

  sqibble_2 <- find_motifs(sq_unt, names_4, c("^VIP", "ONE$", "_"))
  sqibble_2[["found_length"]] <- get_sq_lengths(sqibble_2[["found"]])
  purrr::pwalk(sqibble_2, function(names, sought, found, start, end, found_length) {
    sequence_length <- get_sq_lengths(sq_unt)[which(names == names_4)]
    expect_gte(start, 1)
    expect_gte(end - found_length + 1, 1)
    expect_lte(start + found_length - sequence_length, 1)
    expect_lte(end - sequence_length + 1, 1)
  })

  sqibble_3 <- find_motifs(sq_ami_bsc, names_4, c("XAI", "CX"))
  sqibble_3[["found_length"]] <- get_sq_lengths(sqibble_3[["found"]])
  purrr::pwalk(sqibble_3, function(names, sought, found, start, end, found_length) {
    sequence_length <- get_sq_lengths(sq_ami_bsc)[which(names == names_4)]
    expect_gte(start, 1)
    expect_gte(end - found_length + 1, 1)
    expect_lte(start + found_length - sequence_length, 1)
    expect_lte(end - sequence_length + 1, 1)
  })

  sqibble_4 <- find_motifs(sq_atp, names_3, c("mYmY$", "nsA"))
  sqibble_4[["found_length"]] <- get_sq_lengths(sqibble_4[["found"]])
  purrr::pwalk(sqibble_4, function(names, sought, found, start, end, found_length) {
    sequence_length <- get_sq_lengths(sq_atp)[which(names == names_3)]
    expect_gte(start, 1)
    expect_gte(end - found_length + 1, 1)
    expect_lte(start + found_length - sequence_length, 1)
    expect_lte(end - sequence_length + 1, 1)
  })
})

test_that("index columns can be used to retrieve found subsequence from original sequence", {
  purrr::pwalk(find_motifs(sq_dna_bsc, names_5, "TAG"), function(names, sought, found, start, end) {
    expect_identical(
      bite(sq_dna_bsc[which(names == names_5)], start:end)[[1]],
      found
    )
  })
  purrr::pwalk(find_motifs(sq_ami_bsc, names_4, c("XAI", "CX")), function(names, sought, found, start, end) {
    expect_identical(
      bite(sq_ami_bsc[which(names == names_4)], start:end)[[1]],
      found
    )
  })
  purrr::pwalk(find_motifs(sq_unt, names_4, c("^VIP", "ONE$", "_")), function(names, sought, found, start, end) {
    expect_identical(
      bite(sq_unt[which(names == names_4)], start:end)[[1]],
      found
    )
  })
  purrr::pwalk(find_motifs(sq_atp, names_3, c("mYmY$", "nsA")), function(names, sought, found, start, end) {
    expect_identical(
      bite(sq_atp[which(names == names_3)], start:end)[[1]],
      found
    )
  })
})
