guess_sq_type <- function(x) {
  
}

interpret_type <- function(name) {
  
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


# TODO: sort that rubbish

.guess_sq_type_subtype <- function(x) {
  real_alph <- toupper(.get_real_alph(x))
  .guess_type_subtype_by_alph(real_alph)
}

.guess_type_subtype_by_alph <- function(alph) {
  # TODO: any better idea for it?
  possib_rets <- list(list(type = "dna", is_clean = TRUE),
                      list(type = "rna", is_clean = TRUE),
                      list(type = "ami", is_clean = TRUE),
                      list(type = "dna", is_clean = FALSE),
                      list(type = "rna", is_clean = FALSE),
                      list(type = "ami", is_clean = FALSE))
  for (ret in possib_rets)
    if (all(alph %in% .get_standard_alph(ret[["type"]], ret[["is_clean"]]))) return(ret)
  list(type = "unt", is_clean = NULL)  
}
