#' Find elements which are not suitable for specified type.
#' 
#' @description Finds elements in given sequence not contained in amino acid or
#' nucleotide alphabet.
#' 
#' @template x
#' @template dest_type
#' @template NA_letter
#' @template three-dots
#'
#' @return A list of mismatched elements for every sequence from
#' \code{\link[=sq-class]{sq}} object.
#' 
#' @details
#' Amino acid, DNA and RNA standard alphabets have predefined letters. This
#' function allows the user to check which letters from input sequences are not
#' contained in selected one of these alphabets.
#'
#' Returned list contains a character vector for each input sequence. Each
#' element of a vector is a letter that appear in corresponding sequence and not
#' in the target alphabet.
#'
#' You can check which letters are valid for specified type in
#' \code{\link{alphabet}} documentation.
#'
#' @examples
#' # Creating objects to work on:
#' sq_unt <- sq(c("ACGPOIUATTAGACG","GGATFGHA"), alphabet = "unt")
#' sq_ami <- sq(c("QWERTYUIZXCVBNM","LKJHGFDSAZXCVBN"), alphabet = "ami_ext")
#'
#' # Mismatched elements might be from basic type:
#' find_invalid_letters(sq_ami, "ami_bsc")
#'
#' # But also from type completely unrelated to the current one:
#' find_invalid_letters(sq_unt, "dna_ext")
#' 
#' @family type_functions
#' @seealso \code{\link{alphabet}()}
#' @export
find_invalid_letters <- function(x, dest_type, ...)
  UseMethod("find_invalid_letters")

#' @export
find_invalid_letters.default <- function(x, dest_type, ...)
  stop("method 'find_invalid_letters' isn't implemented for this type of object", call. = FALSE)

#' @rdname find_invalid_letters
#' @export
find_invalid_letters.sq <- function(x, dest_type, ...,
                                    NA_letter = getOption("tidysq_NA_letter")) {
  assert_string(NA_letter, min.chars = 1)
  assert_sq_type(dest_type)
  
  CPP_find_invalid_letters(x, dest_type, NA_letter)
}