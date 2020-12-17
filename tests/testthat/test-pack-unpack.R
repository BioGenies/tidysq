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
sq_unt_6 <- sq(c("POfobyBiY", "", "ZMEepq", "FUasmgh", "iA", "DdG", "voSI",
                 "VLuzQ"),
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
sq_atp_6 <- sq(c("P**fobqoiV", "", "ZMg**pq", "FUbsmgh", "i**", "DDG", "goSI",
                 "VbuG**"),
               alphabet = c(LETTERS, letters, "**"))

local_test_pack_unpack <- function(unpack_format) {
  expect_identical(pack(unpack(sq_dna_bsc, unpack_format), alphabet(sq_dna_bsc)),
                   sq_dna_bsc)
  expect_identical(pack(unpack(sq_dna_ext, unpack_format), alphabet(sq_dna_ext)),
                   sq_dna_ext)
  expect_identical(pack(unpack(sq_rna_bsc, unpack_format), alphabet(sq_rna_bsc)),
                   sq_rna_bsc)
  expect_identical(pack(unpack(sq_rna_ext, unpack_format), alphabet(sq_rna_ext)),
                   sq_rna_ext)
  expect_identical(pack(unpack(sq_ami_bsc, unpack_format), alphabet(sq_ami_bsc)),
                   sq_ami_bsc)
  expect_identical(pack(unpack(sq_ami_ext, unpack_format), alphabet(sq_ami_ext)),
                   sq_ami_ext)

  expect_identical(pack(unpack(sq_unt_2, unpack_format), alphabet(sq_unt_2)),
                   sq_unt_2)
  expect_identical(pack(unpack(sq_unt_4, unpack_format), alphabet(sq_unt_4)),
                   sq_unt_4)
  expect_identical(pack(unpack(sq_unt_6, unpack_format), alphabet(sq_unt_6)),
                   sq_unt_6)

  expect_identical(pack(unpack(sq_atp_2, unpack_format), alphabet(sq_atp_2)),
                   sq_atp_2)
  expect_identical(pack(unpack(sq_atp_3, unpack_format), alphabet(sq_atp_3)),
                   sq_atp_3)
  expect_identical(pack(unpack(sq_atp_4, unpack_format), alphabet(sq_atp_4)),
                   sq_atp_4)
  expect_identical(pack(unpack(sq_atp_5, unpack_format), alphabet(sq_atp_5)),
                   sq_atp_5)
  expect_identical(pack(unpack(sq_atp_6, unpack_format), alphabet(sq_atp_6)),
                   sq_atp_6)
}

# PACK-UNPACK COMPATIBILITY ----
test_that("unpacking and packing returns original sq object for STRING", {
  local_test_pack_unpack("STRING")
})

test_that("unpacking and packing returns original sq object for STRINGS", {
  local_test_pack_unpack("STRINGS")
})

test_that("unpacking and packing returns original sq object for RAWS", {
  local_test_pack_unpack("RAWS")
})

test_that("unpacking and packing returns original sq object for INTS", {
  local_test_pack_unpack("INTS")
})
