#' Find elements which are not suitable for specified type.
#' 
#' Function finds elements in given sequence not matching to amino acid or nucleotide 
#' alphabet. 
#' 
#' @param x a \code{\link{sq}} object to be checked.
#' 
#' @param dest_type a \code{\link{character}} string denoting destination type - it may be 
#' "dna" for \strong{dna} type (DNA), "rna" for \strong{rna} type (RNA) or
#' "ami" for \strong{ami} type (amino acids).  
#'  
#' @return A list of mismatched elements for every sequence from \code{\link{sq}} object.
#' 
#' @details This function allows obtaining list of mismatched elements of sequences from  
#' a \code{\link{sq}} object to amino acid, DNA or RNA alphabet. Output list has number of
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
#' Returned lists for \strong{ami}, \strong{dna} and \strong{rna} \code{sq} objects,
#' if \code{des_type} is specified respectively "ami", "dna" and "rna", will contain
#' only \code{\link[=sq]{NULL}} elements.
#'
#' @seealso \code{\link{sq}} \code{\link{construct_sq}}
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
  assert_sq_type(dest_type)
  
  CPP_find_invalid_letters(x, dest_type, NA_letter)
}