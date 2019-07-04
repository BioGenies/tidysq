#constructing and validating sq

ex_sq1 <- construct_sq("aACccDaaa")
ex_sq2 <- construct_aasq("PmHe")
ex_sq3 <- construct_aasq("QQRdtvIk")
ex_sq4 <- construct_nucsq("ACTTGGAGAGGCT")

tidysq:::validate_sq(ex_sq1)
tidysq:::validate_sq(ex_sq2)
tidysq:::validate_sq(ex_sq3)

tidysq:::validate_ambsq(ex_sq1)
tidysq:::validate_aasq(ex_sq2)
tidysq:::validate_aasq(ex_sq3)

ex_sq5 <- construct_sq("asatgamktr")
tidysq:::validate_ambsq(ex_sq5)
ex_sq6 <- construct_sq("asatgamktr", type = "aa")
tidysq:::validate_aasq(ex_sq6)
ex_sq7 <- construct_sq("asatgamktr", type = "nuc")
tidysq:::validate_nucsq(ex_sq7)

#constructing and validating sqtibble

#at ex_tb1 and ex_tb3 should be warnings
ex_tb1 <- tidysq:::construct_sqtibble(c("name2", "name3", "name4"), list(ex_sq2, ex_sq3, ex_sq4))
ex_tb2 <- tidysq:::construct_sqtibble(c("name2", "name3"), list(ex_sq2, ex_sq3))
ex_tb3 <- tidysq:::construct_sqtibble(c("name1", "name3", "name4"), list(ex_sq1, ex_sq3, ex_sq4))

tidysq:::validate_sqtibble(ex_tb1)
tidysq:::validate_sqtibble(ex_tb2)
tidysq:::validate_sqtibble(ex_tb3)

ex_tb4 <- read_fasta("inst/example_aa.fasta", type = "aa")
tidysq:::validate_sqtibble(ex_tb4)

main_tbaa <- read_fasta("inst/small_example_aa.fasta", type = "aa")
main_tbnuc <- read_fasta("inst/small_example_nuc.fasta", type = "nuc")

clean(main_tbaa)
clean(main_tbaa, only_elements = TRUE)

clean(main_tbnuc)
clean(main_tbnuc, only_elements = TRUE)
