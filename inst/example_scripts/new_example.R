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
