#' Get lengths of sequences in sq object
#' 
#' Function counts number of elements in each sequence in given \code{\link{sq}} object.
#' 
#' @param sq a \code{\link{sq}} object. 
#'  
#' @return A \code{\link{numeric}} vector, where each element gives length of according 
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
  .validate_sq(sq)
  .get_lens(sq)
}