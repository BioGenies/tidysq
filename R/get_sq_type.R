#' Get alphabet of given sq object.
#' 
#' Function returns amino acid, nucleotide or atypical alphabet based on a \code{\link{sq}} 
#' object type. 
#' 
#' @param sq a \code{\link{sq}} object to be recognized. 
#'  
#' @return A character vector of letters of the alphabet.
#' 
#' @details This function allows returning alphabet of \code{sq} object which is a character
#' or numeric vector. The function reads provided \code{\link{sq}} object and determines,
#' which kind of sequences user assigned to a \code{\link{sq}} object (nucleotide, amino acid,
#' atypical or encoded one).
#' 
#' If  \code{\link{sq}} type is \strong{ami} the function returns a set of 20 aminoacids with
#' gap (-) and stop codon (*) letter. If a \code{\link{sq}} contains \strong{nuc} sequences
#' the function returns a set of 5 nucleotides with gap (-) element. If \code{sq} has 
#' additionally \strong{cln} subtype, ambiguous letters are returned as well.
#' 
#' If \code{\link{sq}} type is \strong{unt} or \strong{atp} the function returns a list of 
#' letters present in sequences of a \code{\link{sq}} object.
#' 
#' If type is \strong{enc} a numeric vector of values encoded for letters will be returned
#' (see \code{\link{encode}}).
#' 
#' The details about amino acid and nucleotide 
#' alphabet can be checked in \code{\link{aminoacids_df}} and \code{\link{nucleotides_df}} 
#' respectively. General information about alphabets and types of \code{sq} objects can 
#' be found in \code{\link{sq}} class documentation.
#'
#' @examples 
#' # Creating an object to work on:
#' sq_nucleotides <- construct_sq(c("ACGATTAGACG","GGATA"), type = "nuc")
#' sq_amino_acids <- construct_sq(c("MMVTAAV"), type = "ami")
#' sq_untyped <- construct_sq(c("ACGA&&TTAGACG&"), type = "unt")
#' 
#' # Testing nucleotide sq object with defined type:
#' get_sq_alphabet(sq_nucleotides)
#' 
#' # Testing amino acid sq object with defined type:
#' get_sq_alphabet(sq_amino_acids)
#' 
#' # Testing nucleotide sq object without defined type:
#' get_sq_alphabet(sq_untyped)
#' 
#'   
#' @seealso \code{\link{sq}} \code{\link{construct_sq}} \code{\link{encode}}
#' @export
get_sq_alphabet <- function(sq) {
  validate_sq(sq)
  .get_alph(sq)
}

#' Get type of a sq object
#' 
#' Function checks which type of sequences are contained in \code{\link{sq}} object.
#'  
#' @param sq a \code{\link{sq}} object to be checked. 
#'  
#' @return A \code{\link{character}} string, type of\code{\link{sq}} object - can be one of
#' "ami", "nuc", "unt", "atp" or "enc".
#' 
#' @details This function returns type of sequence from \code{\link{sq}} object.
#' If the type of sequence is \strong{nuc}, \strong{ami}, \strong{unt}, \strong{atp} or 
#' \strong{enc} function returns "nuc", "ami", "unt", "atp" or "enc" respectivetly.
#'  
#' @examples 
#' # Creating an object to work on:
#' sq_nucleotides <- construct_sq(c("ACGATTAGACG","GGATA"), type = "nuc")
#' sq_amino_acids <- construct_sq(c("MMVTAAV"), type = "ami")
#' sq_unt <- construct_sq(c("MMVTAAVXX"), type = "unt")
#' sq_atp <- substitute_letters(sq_amino_acids, c(M = "g1", V = "g2", T = "g1", A = "g3"))
#' sq_encoded <- encode(sq_nucleotides, c(A = 2.3, C = 1.56, T = 0.23, G = 0.28))
#' 
#' # Getting sq type from nucleotide sq object:
#' get_sq_type(sq_nucleotides)
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
get_sq_type <- function(sq) {
  validate_sq(sq)
  .get_sq_type(sq)
}