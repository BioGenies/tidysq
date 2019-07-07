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

####getting sequences types

get_sq_types(sqtbl_1)
get_sq_types(sqtbl_aa)
get_sq_types(sqtbl_nuc)

####dealing with invalid levels
sq_inv1 <- construct_sq("abcdefgh!@#$$%^&&*")
sq_inv2 <- construct_sq("fkdajflka*7^#@")
sqtbl_unt <- construct_sqtibble(c("name1", "name2"), list(sq_inv1, sq_inv2))

get_invalid_levels(sqtbl_unt, "aa")
get_invalid_levels(sqtbl_unt, "aa", only_levels = FALSE)

substitue_invalid_levels(sqtbl_unt, "aa", c(`!` = "A", `@` = "D", `7` = "X", `#` = "F"))

drop_invalid_levels(sqtbl_unt, "aa")

remove_na(sqtbl_1)

get_sq_types(set_sq_types(drop_invalid_levels(sqtbl_unt, "aa"), "aa"))


####cleaning

sqtbl_c1 <- clean(sqtbl_aa)
clean(sqtbl_aa, only_elements = TRUE)
clean(sqtbl_nuc)
clean(sqtbl_nuc, only_elements = TRUE)

####reversing

reverse(sqtbl_aa)
reverse(sqtbl_nuc)

####subsetting sequences

#third example should give warning
bite(sqtbl_aa, 1:3)
bite(sqtbl_nuc, 1:5)
sqtbl_b1 <- bite(sqtbl_4, 1:20)

remove_na(sqtbl_b1)
remove_na(sqtbl_b1, only_elements = TRUE)

####encoding

enc <- c(A = "a", B = "a", C = "a", D = "a", E = "a", F = "b", G = "b", 
         H = "b", I = "c", J = "c", K = "c", L = "c", M = "c", N = "c", 
         O = "c", P = "d", Q = "d", R = "d", S = "d", T = "d", U = "d", 
         V = "d", W = "d", X = "d", Y = "d", Z = "d", `-` = "d")

simplify(sqtbl_2, enc)
get_sq_types(simplify(sqtbl_2, enc))
simplify(sqtbl_4, enc)
simplify(rbind(sqtbl_2, clean(sqtbl_2, only_elements = TRUE)), enc)

#example of reading invalid file, dealing with invalid levels, then cleaning and reducing
## (it was generated with following code):
# paste0(">aa", 1:1000, "\n", sapply(1:1000, function(x)
#   paste0(sample(c(sample(list(c(LETTERS[-24], "#", "-"),
#                               aminoacids_df[!aminoacids_df[["amb"]], "one"]), 1)[[1]],"+"),
#                 sample(4:50, 1), replace = TRUE), collapse = "")), collapse = "\n")

# assume we have fasta file with aa sequences, where - for some unexplainable reason - somebody
# used '#' instead of 'X'; there are also '+' in some sequences, but we don't know, what does it mean

# we want to read file, change '#' to 'X', drop all sequences with '+' (by changing them to NA's
# and then removing them) and reducing alphabet

# here - and only here - we use magrittr's %>% operator for simplicity of code
library(magrittr)

read_fasta("inst/unt_example.fasta") %>%
  substitue_invalid_levels("aa", c(`#` = 'X')) %>%
  drop_invalid_levels("aa") %>%
  remove_na() %>%
  set_sq_types("aa") %>%
  simplify(enc)

####kmers

#this example should give error
count_kmers(sqtbl_2,
            c(1, rep(2, 4), rep(3, 4)),
            list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0), c(1, 1)))
kmers_matrix <- count_kmers(bite(simplify(sqtbl_2, enc), 1:6),
            c(1, rep(2, 4), rep(3, 4)),
            list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0), c(1, 1)))

