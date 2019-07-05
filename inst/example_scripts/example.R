####constructing and validating sq

sq_1 <- construct_sq("aACccDaaa")
sq_2 <- construct_aasq("PmHe")
sq_3 <- construct_aasq("QQRdtvIk")
sq_4 <- construct_nucsq("ACTTGGAGAGGCT")

tidysq:::validate_sq(sq_1)
tidysq:::validate_sq(sq_2)
tidysq:::validate_sq(sq_3)

tidysq:::validate_untsq(sq_1)
tidysq:::validate_aasq(sq_2)
tidysq:::validate_aasq(sq_3)

sq_5 <- construct_sq("asatgamktr")
tidysq:::validate_untsq(sq_5)
sq_6 <- construct_sq("asatgamktr", type = "aa")
tidysq:::validate_aasq(sq_6)
sq_7 <- construct_sq("asatgamktr", type = "nuc")
tidysq:::validate_nucsq(sq_7)

####constructing and validating sqtibble

#at ex_tb1 and ex_tb3 should be warnings
sqtbl_1 <- construct_sqtibble(c("name2", "name3", "name4"), list(sq_2, sq_3, sq_4))
sqtbl_2 <- construct_sqtibble(c("name2", "name3"), list(sq_2, sq_3))
sqtbl_3 <- construct_sqtibble(c("name1", "name3", "name4"), list(sq_1, sq_3, sq_4))

tidysq:::validate_sqtibble(sqtbl_1)
tidysq:::validate_sqtibble(sqtbl_2)
tidysq:::validate_sqtibble(sqtbl_3)

sqtbl_4 <- read_fasta("inst/example_aa.fasta", type = "aa")
tidysq:::validate_sqtibble(sqtbl_4)

####main testing examles

sqtbl_aa <- read_fasta("inst/small_example_aa.fasta", type = "aa")
sqtbl_nuc <- read_fasta("inst/small_example_nuc.fasta", type = "nuc")

####testing simple functions

clean(sqtbl_aa)
clean(sqtbl_aa, only_elements = TRUE)

clean(sqtbl_nuc)
clean(sqtbl_nuc, only_elements = TRUE)

reverse(sqtbl_aa)
reverse(sqtbl_nuc)

#third example should give warning
bite(sqtbl_aa, 1:3)
bite(sqtbl_nuc, 1:5)
bite(sqtbl_4, 1:20)

get_sq_types(sqtbl_1)
get_sq_types(sqtbl_aa)
get_sq_types(sqtbl_nuc)

####encoding
enc <- c(A = "a", B = "a", C = "a", D = "a", E = "a", F = "b", G = "b", 
         H = "b", I = "c", J = "c", K = "c", L = "c", M = "c", N = "c", 
         O = "c", P = "d", Q = "d", R = "d", S = "d", T = "d", U = "d", 
         V = "d", W = "d", X = "d", Y = "d", Z = "d", `-` = "d")

simplify(sqtbl_2, enc)
get_sq_types(simplify(sqtbl_2, enc))
simplify(sqtbl_4, enc)
