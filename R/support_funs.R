#' Get alphabet of given sq object.
#' 
#' Function returns amino acid, nucleotide or atypical alphabet based on a \code{\link{sq}} 
#' object type. 
#' 
#' @param sq a \code{\link{sq}} object to be recognized. 
#'  
#' @return a character vector of letters of the alphabet.
#' 
#' @details This function allows returning alphabet of \code{sq} object which is a character
#' or numeric vector. The function reads provided \code{\link{sq}} object and determines
#' which kind of sequences user assigned to a \code{\link{sq}} object (nucleotide, amino acid,
#' atypical or encoded one).
#' 
#' If  \code{\link{sq}} type is \strong{ami} the function returns a set of 20 aminoacids with
#' gap (-) and stop codon (*) letter. If a \code{\link{sq}} contains \strong{nuc}  sequences
#' the function returns a set of 5 nucleotides with gap (-) element. If \code{sq} has 
#' additionally \strong{cln} subtype, ambiguous letters are returned as well.
#' 
#' If \code{\link{sq}} type is \strong{unt} or \strong{atp} the function returns a list of 
#' letters present in sequences in \code{\link{sq}} object.
#' 
#' If type is \strong{enc} a numeric vector of values encoded for letters will be returned
#' (see \code{\link{encode}}).
#' 
#' The details about amino acid and nucleotide 
#' alphabet can be checked in \code{\link{aminoacids_df}} and \code{\link{nucleotides_df}} 
#' respectively. General informations about alphabets and types of \code{sq} objects can 
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
NULL

#' Get type of a sq object
#' 
#' Function checks which type of sequences are contained in \code{\link{sq}} object.
#'  
#' @param sq a \code{\link{sq}} object to be checked. 
#'  
#' @return a \code{\link{character}} string, type of\code{\link{sq}} object - can be one of
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
NULL


#' Get lengths of sequences in sq object
#' 
#' Function counts number of elements in each sequence in given \code{\link{sq}} object.
#' 
#' @param sq a \code{\link{sq}} object. 
#'  
#' @return a \code{\link{numeric}} vector, where each element gives length of according 
#' sequence from \code{\link{sq}} object.
#' 
#' @details This function allows returning numeric vector of lengths of sequences from
#' \code{\link{sq}} object. The numeric vector is as long as number of sequences present 
#' in \code{\link{sq}} object.
#' The function counts elements in all types of sequences.
#'
#' @examples 
#' # Creating an object to work on:
#' sq_nucleotides <- construct_sq(c("ACGATTAGACG","GGATA"), type = "nuc")
#' sq_amino_acids <- construct_sq(c("MMVTAAV"), type = "ami")
#' 
#' # Counting number of elements in nucleotide sq object with defined type:
#' get_sq_lengths(sq_nucleotides)
#' 
#' # Counting number of elements in amino acid sq object with defined type:
#' get_sq_lengths(sq_amino_acids)
#' 
#' @seealso \code{\link{sq}} \code{\link{construct_sq}}
#' @export
get_sq_lengths <- function(sq) {
  validate_sq(sq)
  .get_lens(sq)
}

#' Test if a sequence is empty
#' 
#' Test a sq object for presence of empty sequences
#' 
#' @param sq \code{\link{sq}} object to be tested
#'  
#' @return a logical vector of the same length as input sq, 
#' indicating which elements are \code{NULL sq}, i.e., an empty sequence
#' of length 0.
#' 
#' @details This function allows identification of empty sequences (that
#' have length 0) represented by the \code{NULL sq} values in the sq object. 
#' It returns a logical for every element of the sq object - \code{TRUE} if
#' its value is \code{NULL sq} and \code{FALSE} otherwise. \code{NULL sq} 
#' values may be introduced as a result of \code{\link{clean}} and 
#' \code{\link{remove_na}} functions. The first one replaces sequences 
#' containing ambiguous elements with \code{NULL sq} values, whereas the 
#' latter replaces \code{NA} values (which may be introduced by 
#' \code{\link{bite}}) with \code{NULL sq}.
#'
#' @examples 
#' # Creating an object to work on:
#' sq <- construct_sq(c("ACGATTAGACG", "", "GACGANTCCAGNTAC"), type = "nuc")
#' 
#' # Testing for presence of empty sequences:
#' is_null_sq(sq)
#' 
#' # Testing for presence of empty sequences after cleaning - sequence 
#' # containing ambiguous elements is replaced by NULL sq:
#' cln_sq <- clean(sq)
#' is_null_sq(cln_sq)
#' 
#' # Testing for presence of empty sequences after using bite
#' # and removing NA.
#' # Extracting letters from first to fifteenth - NA introduced:
#' bitten_sq <- bite(sq, 1:15)
#' # Removing NA:
#' rm_bitten_sq <- remove_na(bitten_sq)
#' # Testing for presence of empty sequences:
#' is_null_sq(rm_bitten_sq)
#' 
#' 
#' @seealso sq clean clnsq
#' @export
is_null_sq <- function(sq) {
  validate_sq(sq)
  unlist(.apply_sq(sq, "char", "none", function(s) length(s) == 0 || 
                     length(s) == 1 && identical(s, "")))
}

#' Find elements which are not suitable for specified type.
#' 
#' Function finds elements in given sequence not matching to amino acid or nucleotide 
#' alphabet. 
#' 
#' @param sq a \code{\link{sq}} object to be checked.
#' 
#' @param dest_type a \code{\link{character}} string denoting destination type - it may be 
#' "nuc" for \strong{nuc} type (nucleotides) or "ami" for \strong{ami} type (amino acids).  
#'  
#' @return a list of mismatched elements for every sequence from \code{\link{sq}} object.
#' 
#' @details This function allows obtaining list of mismatched elements of sequences from  
#' \code{\link{sq}} object to amino acid or nucleotide alphabet. Output list has number of
#' elements equal to length of \code{sq} object and each element is a character vector 
#' of elements that appear in according sequence that does not fit destination type. This 
#' function might be used to find specificly which sequences have which letters - user
#' may want to use this information for example to check input sequeces.
#' 
#' You can check which letters are valid for specified type in \code{\link{sq}} class 
#' documentation.
#' 
#' Mismatched elements might be replaced with other letters or \code{NA} using 
#' \code{\link{substitute_letters}} and then, after removal of all non-standard letters,
#' set type to destinated using \code{\link{typify}}.
#' 
#' Returned lists for \strong{ami} and \strong{nuc} \code{sq} objects, if \code{des_type}
#' is specified respectively "ami" and "nuc", will contain only \code{NULL} elements.
#'
#' @examples
#' # Creating an object to work on:       
#' sq_nucleotides <- construct_sq(c("ACGPOIUATTAGACG","GGATFGHA"))
#' sq_amino_acids <- construct_sq(c("QWERTYUIZXCVBNM","LKJHGFDSAZXCVBN"))
#' 
#' # Creating lists of mismatched elements from nucleotide sq object:
#' get_invalid_letters(sq_nucleotides, "nuc")
#' 
#' # Creating lists of mismatched elements from amino acid sq object:
#' get_invalid_letters(sq_amino_acids, "ami")
#' 
#' @seealso sq construct_sq
#' @export
get_invalid_letters <- function(sq, dest_type) {
  validate_sq(sq)
  .check_type(dest_type, "'dest_type'")
  
  alph <- .get_alph(sq)
  alph_size <- .get_alph_size(alph)
  dest_alph <- .get_standard_alph(dest_type, FALSE)
  dest_alph <- c(dest_alph, tolower(dest_alph))
  
  .apply_sq(sq, "char", "none", function(s) setdiff(s, dest_alph))
}