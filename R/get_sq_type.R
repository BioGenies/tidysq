#' Get alphabet of given sq object.
#' 
#' Function returns amino acid, DNA, RNA or atypical alphabet based on a \code{\link{sq}} 
#' object type. 
#' 
#' @param sq a \code{\link{sq}} object to be recognized. 
#'  
#' @return A character vector of letters of the alphabet.
#' 
#' @details This function allows returning alphabet of \code{sq} object which is a character
#' or numeric vector. The function reads provided \code{\link{sq}} object and determines,
#' which kind of sequences user assigned to a \code{\link{sq}} object (DNA, RNA, amino acid,
#' atypical or encoded one).
#' 
#' If \code{\link{sq}} has \strong{ami} type and \strong{cln} subtype, the function returns
#' the set of 20 aminoacids, the gap (-) and the stop codon (*) letter. If \strong{cln} is
#' missing, ambiguous letters are included as well.
#' 
#' If a \code{\link{sq}} contains \strong{dna} or \strong{rna} sequences and \strong{cln}
#' subtype, the function returns a set of respective 4 nucleotides with gap (-) element.
#' If \strong{cln} is missing, ambiguous letters are included as well.
#' 
#' If \code{\link{sq}} type is \strong{unt} or \strong{atp}, the function returns a list of 
#' letters present in sequences of a \code{\link{sq}} object.
#' 
#' If type is \strong{enc}, a numeric vector of values encoded for letters is returned
#' (see \code{\link{encode}}).
#' 
#' The details about amino acid and nucleotide alphabets can be checked in
#' \code{\link{aminoacids_df}} and \code{\link{nucleotides_df}} respectively. General
#' information about alphabets and types of \code{sq} objects can be found in \code{\link{sq}}
#' class documentation.
#'
#' @examples 
#' # Creating an object to work on:
#' sq_dna <- construct_sq(c("ACGATTAGACG","GGATA"), type = "dna")
#' sq_amino_acids <- construct_sq(c("MMVTAAV"), type = "ami")
#' sq_untyped <- construct_sq(c("ACGA&&TTAGACG&"), type = "unt")
#' 
#' # Testing DNA sq object with defined type:
#' get_sq_alphabet(sq_dna)
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
  .validate_sq(sq)
  .get_alph(sq)
}

#' Get type of a sq object
#' 
#' Function checks which type of sequences are contained in \code{\link{sq}} object.
#'  
#' @param sq a \code{\link{sq}} object to be checked. 
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
#' sq_encoded <- encode(sq_nucleotides, c(A = 2.3, C = 1.56, T = 0.23, G = 0.28))
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
get_sq_type <- function(sq) {
  .validate_sq(sq)
  .get_sq_type(sq)
}