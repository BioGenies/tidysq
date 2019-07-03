ex_sq1 <- construct_sq("aACccDaaa")
ex_sq2 <- construct_aasq("PmHe")
ex_sq3 <- construct_aasq("QQRdtvIk")
ex_sq4 <- construct_nucsq("ACTTGGAGAGGCT")

validate_sq(ex_sq1)
validate_sq(ex_sq2)
validate_sq(ex_sq3)

validate_ambsq(ex_sq1)
validate_aasq(ex_sq2)
validate_aasq(ex_sq3)
