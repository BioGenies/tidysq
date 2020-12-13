# SETUP ----
# 0,1,2,3,4,5,6,7
sq_dna_bsc <- sq(c("TACTGGGCATG", "CAGGTCGGA", "TAGTAGTCCG", "", "ACGGT"),
                 alphabet = "dna_bsc")
sq_dna_ext <- sq(c("NARTYVTCY", "", "ATKCYGDD", "", "DNAKYTD"),
                 alphabet = "dna_ext")
sq_rna_bsc <- sq(c("UAUCAGU-A-GU-CA", "CUG-A-CUGAG-CC", "-CUG-AGAGUA-"),
                 alphabet = "rna_bsc")
sq_rna_ext <- sq(c("", "KBS-UVW-AWWWG", "YGHHH-", "-CRASH", "MND-KUUBV-MY-"),
                 alphabet = "rna_ext")
sq_ami_bsc <- sq(c("ACEH", "PASAI", "MALACCA", "SIAK"),
                 alphabet = "ami_bsc")
sq_ami_ext <- sq(c("OUTLANDISH", "UNSTRUCTURIZED", "FEAR"),
                 alphabet = "ami_ext")

sq_unt_2 <- sq(c("VII", "IVVV", "", "IVVVIVIIVI", "I", "IIVVI", "VIIVIIV",
                 "IVVVVIVVIVIIII"),
               alphabet = "unt")
sq_unt_4 <- sq(c("VIP01", "VIP002", "VIP0004", "MISSING_ONE", "ONEONEONE",
                 "01", "IGNO", "ENGINE_0"),
               alphabet = "unt")

sq_atp_2 <- sq(c("mAnbAnbAmY", "nbAnbAmA", "mYnbA", "", "mAmYmAnbAnbAmAmY",
                 "nbAnbAmAmYnbAmYnbAmYmA", "mAmAmAmYnbA", "mAmAmAmYnbAmY"),
               alphabet = c("mA", "mY", "nbA"))
sq_atp_3 <- sq(c("mAmYmY", "nbAnsAmAmA", "ngYnsAmAnsAnsAngYnbAmY", "nsAmAmAngYmYmYnbA",
                 "mAmYmYnbAmA", "nsA", "nsAnsA", "ngYnsAmYmYnbAnsA"),
               alphabet = c("mA", "mY", "nbA", "nsA", "ngY"))
sq_atp_4 <- sq(c("iaice(?)(?)", "j(?)ajeech(?)cihi", "(?)chi(?)cahi(?)bd",
                 "bef(?)ia", "", "efbihcibfh", "cb(?)(?)iahef", "ag(?)"),
               alphabet = c(letters[1:10], "(?)"))
sq_atp_5 <- sq(c("VPDIN**DVNI**D", "B??PQBOBI**??FNO", "??", "??VNI**??NIV**",
                 "PP??UDO**", "", "PQB??**P", "PABIK**S??IQ??"),
               alphabet = c(LETTERS, "??", "**"))

# PACK-UNPACK COMPATIBILITY ----
test_that("unpacking and packing returns original sq object for STRING", {
  expect_identical(pack(unpack(sq_dna_bsc, "STRING"), alphabet(sq_dna_bsc)),
                   sq_dna_bsc)
  expect_identical(pack(unpack(sq_dna_ext, "STRING"), alphabet(sq_dna_ext)),
                   sq_dna_ext)
  expect_identical(pack(unpack(sq_rna_bsc, "STRING"), alphabet(sq_rna_bsc)),
                   sq_rna_bsc)
  expect_identical(pack(unpack(sq_rna_ext, "STRING"), alphabet(sq_rna_ext)),
                   sq_rna_ext)
  expect_identical(pack(unpack(sq_ami_bsc, "STRING"), alphabet(sq_ami_bsc)),
                   sq_ami_bsc)
  expect_identical(pack(unpack(sq_ami_ext, "STRING"), alphabet(sq_ami_ext)),
                   sq_ami_ext)

  expect_identical(pack(unpack(sq_unt_2, "STRING"), alphabet(sq_unt_2)),
                   sq_unt_2)
  expect_identical(pack(unpack(sq_unt_4, "STRING"), alphabet(sq_unt_4)),
                   sq_unt_4)

  expect_identical(pack(unpack(sq_atp_2, "STRING"), alphabet(sq_atp_2)),
                   sq_atp_2)
  expect_identical(pack(unpack(sq_atp_3, "STRING"), alphabet(sq_atp_3)),
                   sq_atp_3)
  expect_identical(pack(unpack(sq_atp_4, "STRING"), alphabet(sq_atp_4)),
                   sq_atp_4)
  expect_identical(pack(unpack(sq_atp_5, "STRING"), alphabet(sq_atp_5)),
                   sq_atp_5)
})

test_that("unpacking and packing returns original sq object for STRINGS", {
  expect_identical(pack(unpack(sq_dna_bsc, "STRINGS"), alphabet(sq_dna_bsc)),
                   sq_dna_bsc)
  expect_identical(pack(unpack(sq_dna_ext, "STRINGS"), alphabet(sq_dna_ext)),
                   sq_dna_ext)
  expect_identical(pack(unpack(sq_rna_bsc, "STRINGS"), alphabet(sq_rna_bsc)),
                   sq_rna_bsc)
  expect_identical(pack(unpack(sq_rna_ext, "STRINGS"), alphabet(sq_rna_ext)),
                   sq_rna_ext)
  expect_identical(pack(unpack(sq_ami_bsc, "STRINGS"), alphabet(sq_ami_bsc)),
                   sq_ami_bsc)
  expect_identical(pack(unpack(sq_ami_ext, "STRINGS"), alphabet(sq_ami_ext)),
                   sq_ami_ext)

  expect_identical(pack(unpack(sq_unt_2, "STRINGS"), alphabet(sq_unt_2)),
                   sq_unt_2)
  expect_identical(pack(unpack(sq_unt_4, "STRINGS"), alphabet(sq_unt_4)),
                   sq_unt_4)

  expect_identical(pack(unpack(sq_atp_2, "STRINGS"), alphabet(sq_atp_2)),
                   sq_atp_2)
  expect_identical(pack(unpack(sq_atp_3, "STRINGS"), alphabet(sq_atp_3)),
                   sq_atp_3)
  expect_identical(pack(unpack(sq_atp_4, "STRINGS"), alphabet(sq_atp_4)),
                   sq_atp_4)
  expect_identical(pack(unpack(sq_atp_5, "STRINGS"), alphabet(sq_atp_5)),
                   sq_atp_5)
})

test_that("unpacking and packing returns original sq object for RAWS", {
  expect_identical(pack(unpack(sq_dna_bsc, "RAWS"), alphabet(sq_dna_bsc)),
                   sq_dna_bsc)
  expect_identical(pack(unpack(sq_dna_ext, "RAWS"), alphabet(sq_dna_ext)),
                   sq_dna_ext)
  expect_identical(pack(unpack(sq_rna_bsc, "RAWS"), alphabet(sq_rna_bsc)),
                   sq_rna_bsc)
  expect_identical(pack(unpack(sq_rna_ext, "RAWS"), alphabet(sq_rna_ext)),
                   sq_rna_ext)
  expect_identical(pack(unpack(sq_ami_bsc, "RAWS"), alphabet(sq_ami_bsc)),
                   sq_ami_bsc)
  expect_identical(pack(unpack(sq_ami_ext, "RAWS"), alphabet(sq_ami_ext)),
                   sq_ami_ext)

  expect_identical(pack(unpack(sq_unt_2, "RAWS"), alphabet(sq_unt_2)),
                   sq_unt_2)
  expect_identical(pack(unpack(sq_unt_4, "RAWS"), alphabet(sq_unt_4)),
                   sq_unt_4)

  expect_identical(pack(unpack(sq_atp_2, "RAWS"), alphabet(sq_atp_2)),
                   sq_atp_2)
  expect_identical(pack(unpack(sq_atp_3, "RAWS"), alphabet(sq_atp_3)),
                   sq_atp_3)
  expect_identical(pack(unpack(sq_atp_4, "RAWS"), alphabet(sq_atp_4)),
                   sq_atp_4)
  expect_identical(pack(unpack(sq_atp_5, "RAWS"), alphabet(sq_atp_5)),
                   sq_atp_5)
})

test_that("unpacking and packing returns original sq object for INTS", {
  expect_identical(pack(unpack(sq_dna_bsc, "INTS"), alphabet(sq_dna_bsc)),
                   sq_dna_bsc)
  expect_identical(pack(unpack(sq_dna_ext, "INTS"), alphabet(sq_dna_ext)),
                   sq_dna_ext)
  expect_identical(pack(unpack(sq_rna_bsc, "INTS"), alphabet(sq_rna_bsc)),
                   sq_rna_bsc)
  expect_identical(pack(unpack(sq_rna_ext, "INTS"), alphabet(sq_rna_ext)),
                   sq_rna_ext)
  expect_identical(pack(unpack(sq_ami_bsc, "INTS"), alphabet(sq_ami_bsc)),
                   sq_ami_bsc)
  expect_identical(pack(unpack(sq_ami_ext, "INTS"), alphabet(sq_ami_ext)),
                   sq_ami_ext)

  expect_identical(pack(unpack(sq_unt_2, "INTS"), alphabet(sq_unt_2)),
                   sq_unt_2)
  expect_identical(pack(unpack(sq_unt_4, "INTS"), alphabet(sq_unt_4)),
                   sq_unt_4)

  expect_identical(pack(unpack(sq_atp_2, "INTS"), alphabet(sq_atp_2)),
                   sq_atp_2)
  expect_identical(pack(unpack(sq_atp_3, "INTS"), alphabet(sq_atp_3)),
                   sq_atp_3)
  expect_identical(pack(unpack(sq_atp_4, "INTS"), alphabet(sq_atp_4)),
                   sq_atp_4)
  expect_identical(pack(unpack(sq_atp_5, "INTS"), alphabet(sq_atp_5)),
                   sq_atp_5)
})
