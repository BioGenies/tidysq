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
#' @return A list of mismatched elements for every sequence from \code{\link{sq}} object.
#' 
#' @details This function allows obtaining list of mismatched elements of sequences from  
#' a \code{\link{sq}} object to amino acid or nucleotide alphabet. Output list has number of
#' elements equal to length of \code{sq} object and each element is a character vector 
#' of elements that appear in according sequence that does not fit destination type. This 
#' function might be used to find specifically, which sequences have letters - user
#' may want to use this information for example to check input sequences.
#' 
#' You can check, which letters are valid for specified type in \code{\link{sq}} class 
#' documentation.
#' 
#' Mismatched elements might be replaced with other letters or \code{NA} using 
#' \code{\link{substitute_letters}} and then, after removal of all non-standard letters,
#' set type to destinated using \code{\link{typify}}.
#' 
#' Returned lists for \strong{ami} and \strong{nuc} \code{sq} objects, if \code{des_type}
#' is specified respectively "ami" and "nuc", will contain only \code{\link[=sq]{NULL}} elements.
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
#' @seealso \code{\link{sq}} \code{\link{construct_sq}}
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