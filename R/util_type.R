obtain_alphabet <- function(x, sample_size = 4096, 
                            NA_letter = getOption("tidysq_NA_letter")) {
  CPP_obtain_alphabet(x, sample_size, NA_letter)
}

guess_standard_alphabet <- function(alph,
                                    NA_letter = getOption("tidysq_NA_letter")) {
  CPP_guess_standard_alph(alph, NA_letter)
}

interpret_type <- function(name) {
  # TODO: improve; just improve
  switch(name, 
         dna_bsc = "dna_bsc",
         dna_ext = "dna_ext",
         rna_bsc = "rna_bsc",
         rna_ext = "rna_ext",
         ami_bsc = "ami_bsc",
         ami_ext = "ami_ext",
         unt = "unt")
}

type_as_class <- function(type)
  paste0("sq_", type)

#' Get type of a sq object
#' 
#' Function checks which type of sequences are contained in \code{\link{sq}} object.
#'  
#' @param x a \code{\link{sq}} object to be checked.
#'  
#' @return A \code{\link{character}} string, type of\code{\link{sq}} object - can be one of
#' "ami", "dna", "rna", "unt", "atp" or "enc".
#' 
#' @details This function returns type of sequence from \code{\link{sq}} object.
#' If the type of sequence is \strong{dna}, \strong{rna}, \strong{ami}, \strong{unt},
#' \strong{atp} or \strong{enc} function returns "dna", "rna", "ami", "unt", "atp" or
#' "enc" respectivetly.
#'  
#' @examples 
#' # Creating an object to work on:
#' sq_dna <- construct_sq(c("ACGATTAGACG","GGATA"), type = "dna")
#' sq_rna <- construct_sq(c("CCUUACGGC","UUCAAAGCU"), type = "rna")
#' sq_amino_acids <- construct_sq(c("MMVTAAV"), type = "ami")
#' sq_unt <- construct_sq(c("MMVTAAVXX"), type = "unt")
#' sq_atp <- substitute_letters(sq_amino_acids, c(M = "g1", V = "g2", T = "g1", A = "g3"))
#' sq_encoded <- encode(sq_dna, c(A = 2.3, C = 1.56, T = 0.23, G = 0.28))
#' 
#' # Getting sq type from DNA sq object:
#' get_sq_type(sq_dna)
#' 
#' # Getting sq type from RNA sq object:
#' get_sq_type(sq_rna)
#' 
#' # Getting sq type from amino acid  sq object:
#' get_sq_type(sq_amino_acids)
#' 
#' # Getting sq type from untyped sq object:
#' get_sq_type(sq_unt)
#' 
#' # Getting sq type from atypical sq object:
#' get_sq_type(sq_atp)
#' 
#' # Getting sq type from encoded sq object:
#' get_sq_type(sq_encoded)
#' 
#' @seealso \code{\link{sq}} \code{\link{construct_sq}}
#' @export
get_sq_type <- function(x) {
  # TODO: a generic, maybe?
  assert_class(x, "sq")
  vec_ptype_abbr(x)
}
