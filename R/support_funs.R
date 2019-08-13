#' Returns full alphabet of type of given sequence.
#' 
#' Function returns amino acid or nucleotide alphabet based on sequence type. 
#' 
#' @param sq \code{\link{sq}} object to be recognized. 
#'  
#' @return a character vector of amino acid or nucleotide alphabet elements.
#' 
#' @details This function allows returning character vector of amino acid or 
#' nucleotide alphabet. The function read provided \code{\link{sq}} object and determines
#' which kind of sequence user assigned to a \code{\link{sq}} object (nucleotide or amino acid).
#' If  \code{\link{sq}} contains an amino acid sequence the function returns set of 20 aminoacids with
#' gap (-) and any amino acid (*) element. If a \code{\link{sq}} contains
#' a nucleotide sequence the function returns set of 5 nucleotides with gap (-) element.
#' If \code{\link{sq}} type is defined as "unt" function returns list of elements present
#' in \code{\link{sq}} object. The details about amino acid and nucleotide alphabet can be
#' checked in \code{\link{aminoacids_df}} and \code{\link{nucleotides_df}} respectively.
#'
#' @examples 
#' # Creating an object to work on:
#' sq_nucleotides <- construct_sq(c("ACGATTAGACG","GGATA"), type = "nuc")
#' sq_amino_acids <- construct_sq(c("MMVTAAV"), type = "ami")
#' sq_without_type <- construct_sq(c("ACGATTAGACG"))
#' 
#' # Testing nucleotide sq object with defined type:
#' get_sq_alphabet(sq_nucleotides)
#' 
#' # Testing amino acid sq object with defined type:
#' get_sq_alphabet(sq_amino_acids)
#' 
#' # Testing nucleotide sq object without defined type:
#' get_sq_alphabet(sq_without_type)
#' 
#'   
#' @seealso sq construct_sq
#' @export
get_sq_alphabet <- function(sq) {
  validate_sq(sq)
  .get_alph(sq)
}
NULL

#' Checks type of the sequence 
#' 
#' Function is checking which type of sequence contains \code{\link{sq}} object.
#'  
#' @param sq \code{\link{sq}} object to be checked. 
#'  
#' @return a type string and its name defined in \code{\link{sq}} object.
#' 
#' @details This function returns type of sequence from \code{\link{sq}} object.
#' If the type of sequence is "nuc", "ami" or "unt"  function returns 
#' (nucsq) "nuc", (amisq) "ami" or (untsq) "unt" respectivetly.
#'  
#' @examples 
#' # Creating an object to work on:
#' sq_nucleotides <- construct_sq(c("ACGATTAGACG","GGATA"), type = "nuc")
#' sq_amino_acids <- construct_sq(c("MMVTAAV"), type = "ami")
#' sq_unt <- construct_sq(c("MMVTAAVXX"), type = "unt")
#' 
#' # Gettig sequence type from nucleotide sq object:
#' get_sq_type(sq_nucleotides)
#' 
#' # Gettig sequence type from amino acid  sq object:
#' get_sq_type(sq_amino_acids)
#' 
#' # Gettig sequence type from unt sq object:
#' get_sq_type(sq_unt)
#' 
#' 
#' @seealso sq construct_sq
#' @export
get_sq_type <- function(sq) {
  validate_sq(sq)
  .get_sq_type(sq)
}
NULL


#' Checks how long the given sequence is.
#' 
#' Function is counting number of elements from given sequence. 
#' 
#' @param sq \code{\link{sq}} object to be counted. 
#'  
#' @return a numeric vector of counted elements for every sequence from \code{\link{sq}} object.
#' 
#' @details This function allows returning numeric vector of counted elements from
#' \code{\link{sq}} object. Function counts elements for every sequence from \code{\link{sq}} object.
#' The numeric vector is long as number of sequences present in \code{\link{sq}} object.
#' The function counts elements in nucleotide as well as amino acid sequences.
#'
#' @examples 
#' # Creating an object to work on:
#' sq_nucleotides <- construct_sq(c("ACGATTAGACG","GGATA"), type = "nuc")
#' sq_amino_acids <- construct_sq(c("MMVTAAV"), type = "ami")
#' sq_unt <- construct_sq(c("MMVTAAVXX"), type = "unt")
#' sq_without_type <- construct_sq(c("ACGATTAGACG"))
#' 
#' # Counting number of elements in nucleotide sq object with defined type:
#' get_sq_lenghts(sq_nucleotides)
#' 
#' # Counting number of elements in amino acid sq object with defined type:
#' get_sq_lenghts(sq_amino_acids)
#' 
#' # Counting number of elements in unt sq object with defined type:
#' get_sq_lenghts(sq_unt)
#' 
#' # Counting number of elements in nucleotide sq object without defined type:
#' get_sq_lenghtst(sq_without_type)
#' 
#' 
#' @seealso sq construct_sq 
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

#' Finds elements which are not present in nucleotide or amino acid alphabet.
#' 
#' Function finds elements in given sequence not matching to amino acid or nucleotide 
#' alphabet. 
#' 
#' @param sq \code{\link{sq}} object to be checked.
#' 
#' @param dest_type indices specifying alphabet. Its may be 
#' \code{nuc} for nucleotide alphabet or \code{ami} for amino acid alphabet.  
#'  
#' @return a list of mismatched elements for every sequence from \code{\link{sq}} object.
#' 
#' @details This function allows returning list of mismatched elements from  
#' \code{\link{sq}} object to amino acid or nucleotide alphabet. Function creates list
#' for each sequence from \code{\link{sq}} object and append to each list  element not 
#' matching to any alphabet. 
#'
#' @examples
#' # Creating an object to work on:       
#' sq_nucleotides <- construct_sq(c("ACGPOIUATTAGACG","GGATFGHA"), type = "nuc")
#' sq_amino_acids <- construct_sq(c("QWERTYUIZXCVBNM","LKJHGFDSAZXCVBN"), type = "ami")
#' 
#' # Creating lists of mismatched elements from nucleotide sq object:
#' get_invalid_letters(sq_nucleotides)
#' 
#' # Creating lists of mismatched elements from amino acid sq object:
#' get_invalid_letters(sq_amino_acids)
#' 
#' @seealso sq construct_sq
#' @export
get_invalid_letters <- function(sq, dest_type) {
  validate_sq(sq)
  if (missing(dest_type) || 
     !(dest_type %in% c("ami", "nuc"))) {
    stop("'dest_type' should be either 'ami' or 'nuc'")
  }
  
  alph <- .get_alph(sq)
  alph_size <- .get_alph_size(alph)
  na_char <- .get_na_char()
  dest_alph <- .get_standard_alph(dest_type, FALSE)
  dest_alph <- c(dest_alph, tolower(dest_alph))
  
  .apply_sq(sq, "char", "none", function(s) setdiff(s, dest_alph))
}


