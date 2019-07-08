#####################################################################
#constructing sq:

sq_1 <- construct_sq("PmHe", type = "aa") #subclass : aasq (type : aa)
sq_2 <- construct_sq("QQRdtvIk", type = "aa")
sq_3 <- construct_sq("ACTTGGAGAGGCT", type = "nuc") #subclass : nucsq (type : nuc)

construct_sq("##&*^#!r2q431948", type = "nuc") # there will be error

sq_4 <- construct_sq("aACccDaaa") #subclass : untsq (class : unt)

sq_5 <- construct_sq(c("a", "B", "c", "D", "e", "F"))

#TO DO: constructing aasq from three-letters aminoacids names

#printing:
sq_1
sq_2
sq_3
sq_4

#validating:
#QUESTION: should this method be exported?
tidysq:::validate_sq(sq_1)
tidysq:::validate_sq(sq_2)
tidysq:::validate_sq(sq_3)
tidysq:::validate_sq(sq_4)

tidysq:::validate_aasq(sq_2)

#constructing sqtibble:
sqtbl_1 <- construct_sqtibble(c("name1", "name2", "name3"), list(sq_1, sq_2, sq_3)) #warning
sqtbl_2 <- construct_sqtibble(c("name1", "name2"), list(sq_1, sq_2))
sqtbl_3 <- construct_sqtibble(c("name1", "name3", "name4"), list(sq_1, sq_3, sq_4)) #warnings

sqtbl_1[["sq"]]

#QUESTION: should this be exported?
tidysq:::validate_sqtibble(sqtbl_1)


#TO DO: constructing with given string of named vector of strings, like 
# c(name1 = "AACVVDA", name2 ="ASDASDA")
# >> QUESTION: should we allow constructing sq objects from vetors like c("A", "B", "C")?

#TO DO: not allowing keeping different types of sq objects in sqtibble
# >> QUESTION: for sure? 

#QUESTION: is printing ok?
#TO DO: modify printing -add information about type of sq

#####################################################################
#exported data:
data("aaprop")
data("aminoacids_df")
data("nucleotides_df")

#QUESTION - is it correct?

#####################################################################
#reading files:

sqtbl_long <- read_fasta("inst/example_aa.fasta", type = "aa")
sqtbl_aa <- read_fasta("inst/small_example_aa.fasta", type = "aa")
sqtbl_nuc <- read_fasta("inst/small_example_nuc.fasta", type = "nuc")
#QUESTION: what about comments in files?

#getting types:
get_sq_types(sqtbl_1)
get_sq_types(sqtbl_aa)
get_sq_types(sqtbl_nuc)
# NOTE: if we remove different types of sq, then it won't be useful

#####################################################################
#dealing with improper data:
sq_inv1 <- construct_sq("abcdefgh!@#$$%^&&*")
sq_inv2 <- construct_sq("fkdajflka*7^#@")
sqtbl_unt <- construct_sqtibble(c("name1", "name2"), list(sq_inv1, sq_inv2))

get_invalid_levels(sqtbl_unt, "aa")
get_invalid_levels(sqtbl_unt, "aa", only_levels = FALSE)

sqtbl_unt
substitue_invalid_levels(sqtbl_unt, "aa", c(`!` = "A", `@` = "D", `7` = "X", `#` = "F"))

drop_invalid_levels(sqtbl_unt, "aa")
drop_invalid_levels(sqtbl_unt, "aa", replacement = "X")

remove_na(drop_invalid_levels(sqtbl_unt, "aa"))
remove_na(drop_invalid_levels(sqtbl_unt, "aa"), only_elements = TRUE)
#QUESTION: what about the name? should this be just an alias for overloaded na.omit()?

set_sq_types(remove_na(drop_invalid_levels(sqtbl_unt, "aa"), only_elements = TRUE), "aa")

#example:
library(magrittr)

read_fasta("inst/unt_example.fasta") %>%
  substitue_invalid_levels("aa", c(`#` = 'X', `+`="S")) %>%
  drop_invalid_levels("aa") %>%
  remove_na() %>%
  set_sq_types("aa")

#####################################################################
#cleaning:
clean(sqtbl_aa)
clean(sqtbl_aa, only_elements = TRUE)
clean(sqtbl_nuc)
clean(sqtbl_nuc, only_elements = TRUE)

clean(sqtbl_nuc)[["sq"]] #keeps info about being cleaned
class(clean(sqtbl_nuc)[["sq"]][[1]]) #another subclass added
class(clean(sqtbl_aa)[["sq"]][[1]]) 
#reversing:
reverse(sqtbl_aa)
reverse(sqtbl_nuc)

#biting:
bite(sqtbl_aa, 1:3)
bite(sqtbl_aa, -2)
bite(sqtbl_nuc, 1:5)
bite(sqtbl_long, 1:20)
#QUESTION: what about the name?

#####################################################################
enc <- c(A = "a", B = "a", C = "a", D = "a", E = "a", F = "b", G = "b", 
         H = "b", I = "c", J = "c", K = "c", L = "c", M = "c", N = "c", 
         O = "c", P = "d", Q = "d", R = "d", S = "d", T = "d", U = "d", 
         V = "d", W = "d", X = "d", Y = "d", Z = "d", `-` = "d")

simplify(sqtbl_aa, enc) #new sq subclass: sqsim (type: sim)
get_sq_types(simplify(sqtbl_aa, enc))
simplify(sqtbl_long, enc)
#warning expected due to cleaned and uncleaned sequences at once; 
# >> if we won't allow it, it won't be necessary
simplify(rbind(sqtbl_2, clean(sqtbl_2, only_elements = TRUE)), enc) 

#QUESTION: what about the names? simplify for reducing alphabet, encoding for getting 
#encoded sequences?

#####################################################################
count_kmers(sqtbl_aa,
            c(1, rep(2, 4), rep(3, 4)),
            list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0), c(1, 1)))
kmers_matrix <- count_kmers(bite(simplify(sqtbl_aa, enc), 1:6),
                            c(1, rep(2, 4), rep(3, 4)),
                            list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0), c(1, 1)))
