library(dplyr)

###### ----- creating and reading data ----- ######
###################################################

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

tidysq:::validate_sq(sq_1, type = "nuc")
tidysq:::validate_sq(sq_2, type = "nuc")
tidysq:::validate_sq(sq_3, type = "ami")
tidysq:::validate_sq(sq_4, type = "ami")
tidysq:::validate_sq(sq_5, type = "unt")

### reading and writing fasta

sqtbl_ami <- read_fasta("inst/small_example_aa.fasta", type = "ami")
sqtbl_nuc <- read_fasta("inst/small_example_nuc.fasta", type = "nuc")
sqtbl_long <- read_fasta("inst/example_aa.fasta", type = "ami")
# TODO: FEATURE - guessing of sequences type

write_fasta(sq = sqtbl_long[["sq"]], 
            name = sqtbl_long[["name"]], 
            file = "inst/save.fasta")
# TODO: INPUT VALIDATION - not allowing saving sequences with NA's

###### ----- dealing with data ----- ######
###########################################

### get_invalid_letters

get_invalid_letters(sq_5, "ami") 

### substitute_letters

substitute_letters(sq_4, c(F = "D", A = ";", D = ";"))
substitute_letters(sq_4, c(F = "D", A = NA, D = NA))
sqtbl_ami %>% mutate(subs_1 = substitute_letters(sq, c(G = "L", V = ";")),
                     subs_2 = substitute_letters(sq, c(P = "G", K = NA, X = NA, V = NA)))
#QUESTION: which character would be best for printing NA's?

### remove_na function
remove_na(sq_4)
remove_na(sq_4, only_elements = TRUE)

### typify

typify(substitute_letters(sq_5, 
                          c(`2` = "A", `4` = "B", `3` = "X",`;` = "X", `'` = "X", `9` = "X")), 
       "ami")

tibble(sq = sq_5[-1]) %>% mutate(subst = substitute_letters(sq, c(`2` = "A", `4` = "B", `3` = NA,
                                                    `;` = NA, `'` = NA, `9` = NA)),
                   removed = remove_na(subst),
                   typed = typify(removed, "ami"))

###### ----- operating on data ----- ######
###########################################

### clean function

sqtbl_ami %>% mutate(cleaned = clean(sq))
sqtbl_ami %>% mutate(cleaned = clean(sq, only_elements = TRUE))
sqtbl_nuc %>% mutate(cleaned = clean(sq))
sqtbl_nuc %>% mutate(cleaned = clean(sq, only_elements = TRUE))

### reverse function

reverse(sq_4)
sqtbl_ami %>% mutate(reversed = reverse(sq))
sqtbl_nuc %>% mutate(reversed = reverse(sq))

### bite function

bite(sq_4, 1:4)
sqtbl_ami %>% mutate(bitten = bite(sq, 1:3))
sqtbl_ami %>% mutate(bitten = bite(sq, -5:-8))
sqtbl_ami %>% mutate(bitten = bite(sq, 1:20))

### remove_na comes back

sqtbl_ami %>% mutate(bitten = remove_na(bite(sq, 1:20), only_elements = TRUE))
sqtbl_ami %>% mutate(bitten = bite(sq, 1:15),
                     na_removed = remove_na(bitten),
                     reveresed = reverse(na_removed))

### complement
complement(clean(sqtbl_nuc %>% pull(sq)))  #don't need to specify if is_dna
complement(clean(sqtbl_nuc %>% pull(sq)), is_dna = TRUE) #it is optional in this case
sqtbl_nuc %>% mutate(rnaed = substitute_letters(sq, c(T = "U")),
                     nuced = typify(rnaed, "nuc"),
                     cleaned = clean(nuced, only_elements = TRUE),
                     compl_1 = complement(cleaned),
                     compl_2 = complement(cleaned, is_dna = FALSE)) #as well here
tibble(sq = construct_sq(c("TGCGCGT", "TGC", "CTG"), type = "nuc")) %>%
  mutate(sq_2 = complement(clean(sq))) #here, as there is no 'A' in sequences, you don't need specification

### simplfy
enc <- c(A = "a", B = "a", C = "a", D = "a", E = "a", F = "b", G = "b", 
         H = "b", I = "c", J = "c", K = "c", L = "c", M = "c", N = "c", 
         O = "c", P = "d", Q = "d", R = "d", S = "d", T = "d", U = "d", 
         V = "d", W = "d", X = "d", Y = "d", Z = "d", `-` = "d", `*` = "d")

simplify(sqtbl_ami %>% pull(sq), enc)
sqtbl_ami %>% mutate(simpl = simplify(sq, enc))
#TODO: FEATURE - find a better way to write encodings

### encode
data("aaprop")
enc2 <- aaprop[1,]

sqtbl_ami[["sq"]] %>% clean(only_elements = TRUE) %>% encode(enc2) %>% as.matrix
sqtbl_ami %>% mutate(enc = encode(clean(sq, only_elements = TRUE), enc2))
sqtbl_long %>% mutate(enc = encode(clean(sq, only_elements = TRUE), enc2))

###### ----- overloaded methods for sq class ----- ######
#########################################################

### as.character

as.character(sq_1)
as.character(sq_2)
as.character(sq_5)
as.character(bite(sqtbl_nuc[["sq"]], 1:20))

### is.sq

is.sq(sq_1)
is.sq(sq_2)
is.sq(sq_3)
is.sq(sq_4)
is.sq(sq_5)

is.amisq(sq_1)
is.amisq(sq_2)
is.amisq(sq_3)
is.amisq(sq_4)
is.amisq(sq_5)

is.nucsq(sq_1)
is.nucsq(sq_2)
is.nucsq(sq_3)
is.nucsq(sq_4)
is.nucsq(sq_5)

is.untsq(sq_1)
is.untsq(sq_2)
is.untsq(sq_3)
is.untsq(sq_4)
is.untsq(sq_5)

is.simsq(sqtbl_ami %>% pull(sq))
is.simsq(simplify(sqtbl_ami %>% pull(sq), enc))

is.atpsq(sq_1)
is.atpsq(substitute_letters(sq_1, c(A = 'F')))

### ==

sq_1 == "ACTAGAGTGATAGTAGGAGTAGA"
sq_2 == c("CCCCC","TTT", "ACTG", "ACTC")
sq_3 == "AajsfdjLKFAJkajd"
sq_4 == c("fafasfasfFSA", "ygagayagfa", "adsDaf")
sq_5 == c("afsfd", "q243faadfa", "afsw34gesfv", "adfq2", "fasfas", "g'qp9u2r3'b;")

# for 'ami' and 'nuc' sequences size doesn't matter, for others - does
sq_1 == tolower("ACTAGAGTGATAGTAGGAGTAGA")
sq_4 == tolower(c("fafasfasfFSA", "ygagayagfa", "adsDaf"))
sq_5 == toupper(c("afsfd", "q243faadfa", "afsw34gesfv", "adfq2", "fasfas", "g'qp9u2r3'b;"))

sq_1 == sq_1
sq_1 == sq_2
sqtbl_ami[["sq"]][-3] == clean(sqtbl_ami[["sq"]])[-3]

###### ----- filtering functions ----- ######
#############################################

### operator %has%

# for ami and nuc sq it translates some letters accordingly to standard; it also treats all like uppers

sq_ami <- (sqtbl_ami %>% pull("sq"))

sq_ami %has% "GG"
sq_ami %has% "n"
sq_ami %has% "J" # translates J into L, I or J

sq_ami %has% c("K", "P", "Q")
sq_ami %has% "K" ||
  sq_ami %has% "P" ||
  sq_ami %has% "Q"
sq_ami %has% "KPQ"

sq_ami %has% "IVYKpvdLSKVT"

(sqtbl_nuc %>% pull("sq")) %has% "GG"
(sqtbl_nuc %>% pull("sq")) %has% "GtaTGCT"
(sqtbl_nuc %>% pull("sq")) %has% "CN" # translates N into any aminoacid
(sqtbl_nuc %>% pull("sq")) %has% c("GC", "at")

construct_sq(c("CTGA-N", "ACTGH", "SD"), type = "nuc") %has% "AN"
construct_sq(c("CTGA-N", "ACTGH", "SD"), type = "nuc") %has% "A-" # N is any but gap

sq_5 %has% "faa"
sq_5 %has% "af"
sq_5 %has% c("a", "2")
sq_5 %has% c("^a", "s")

(sqtbl_long %>% pull("sq") %>% simplify(enc)) %has% "acda"

sqtbl_long %>%
  filter(sq %has% c("KLV", "^S", "HxxxxxF"))

sqtbl_long %>%
  filter(sq %has% c("^D", "A$"))

### is_null_sq

is_null_sq(clean(sq_ami))


### find_motifs

find_motifs(sqtbl_long[["sq"]], sqtbl_long[["name"]], c("AS"))
find_motifs(sqtbl_long[["sq"]], sqtbl_long[["name"]], c("X", "DF"))
find_motifs(sqtbl_long[["sq"]], sqtbl_long[["name"]], c("XXX"))
find_motifs(sqtbl_long[["sq"]], sqtbl_long[["name"]], c("^D"))
#QUESTION: is name required?

### more advanced example:

read_fasta("inst/unt_example.fasta") %>%
  mutate(sq = sq %>% 
           substitute_letters(c(`#` = "X", `+` = NA)) %>% 
           remove_na()) %>%
  mutate(sq = sq %>% 
           typify("ami") %>%
           clean()) %>%
  filter(!is_null_sq(sq)) %>%
  filter(lengths(sq) > 6) %>%
  mutate(sq = sq %>% 
           simplify(enc) %>%
           bite(1:18) %>%
           remove_na(only_elements = TRUE)) %>%
  filter(sq %has% "cccc")

