library(dplyr)

### sq object - creating, printing and validation

sq_1 <- construct_sq("ACTAGAGTGATAGTAGGAGTAGA", type = "nuc")
sq_2 <- construct_sq(c("CCCCC","TTT", "ACTG", "ACTC"), type = "nuc")
sq_3 <- construct_sq("AajsfdjLKFAJkajd", type = "ami")
sq_4 <- construct_sq(c("fafasfasfFSA", "ygagayagfa", "adsDaf"), type = "ami")
sq_5 <- construct_sq(c("afsfd", "q243faadfa", "afsw34gesfv", "adfq2", "fasfas", "g'qp9u2r3'b;"))

sq_1
sq_2
sq_3
sq_4
sq_5

tidysq:::validate_sq(sq_1)
tidysq:::validate_sq(sq_2)
tidysq:::validate_sq(sq_3)
tidysq:::validate_sq(sq_4)
tidysq:::validate_sq(sq_5)

tidysq:::validate_nucsq(sq_1)
tidysq:::validate_nucsq(sq_2)
tidysq:::validate_amisq(sq_3)
tidysq:::validate_amisq(sq_4)
tidysq:::validate_untsq(sq_5)

### sqtibble object - creating, printing and validation

sqtbl_1 <- construct_sqtibble(sq_1)
sqtbl_2 <- construct_sqtibble(sq_2)
sqtbl_3 <- construct_sqtibble(sq_3)
sqtbl_4 <- construct_sqtibble(sq_4)
sqtbl_5 <- construct_sqtibble(sq_5)
sqtbl_6 <- construct_sqtibble(c("asfawfaf", "kifachbjhj", "jasfbfka"), type = "ami")
sqtbl_7 <- construct_sqtibble(c(nuc1 = "CCTGAGTA", nuc2 = "TAGTCTAGTAGA"), type = "nuc")
sqtbl_8 <- construct_sqtibble(c("akuu", "2q43hor", "A981O  hkJA"), name = c("n1", "n2", "n3"))

sqtbl_1
sqtbl_2
sqtbl_3
sqtbl_4
sqtbl_5
sqtbl_6
sqtbl_7
sqtbl_8

tidysq:::validate_sqtibble(sqtbl_1)
tidysq:::validate_sqtibble(sqtbl_2)
tidysq:::validate_sqtibble(sqtbl_3)
tidysq:::validate_sqtibble(sqtbl_4)
tidysq:::validate_sqtibble(sqtbl_5)
tidysq:::validate_sqtibble(sqtbl_6)
tidysq:::validate_sqtibble(sqtbl_7)
tidysq:::validate_sqtibble(sqtbl_8)

### reading fasta

sqtbl_ami <- read_fasta("inst/small_example_aa.fasta", type = "ami")
sqtbl_nuc <- read_fasta("inst/small_example_nuc.fasta", type = "nuc")
sqtbl_long <- read_fasta("inst/example_aa.fasta", type = "ami")

tidysq:::validate_sqtibble(sqtbl_ami)
tidysq:::validate_sqtibble(sqtbl_nuc)
tidysq:::validate_sqtibble(sqtbl_long)

### clean function

sqtbl_ami %>% mutate(cleaned = clean(sq))
sqtbl_ami %>% mutate(cleaned = clean(sq, only_elements = TRUE))
sqtbl_nuc %>% mutate(cleaned = clean(sq))
sqtbl_nuc %>% mutate(cleaned = clean(sq, only_elements = TRUE))

### reverse function

sqtbl_ami %>% mutate(reversed = reverse(sq))
sqtbl_nuc %>% mutate(reversed = reverse(sq))

