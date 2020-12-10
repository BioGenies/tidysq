#' Reverse sequence
#' 
#' @description Reverse given list of sequences.
#' 
#' @template x
#' @template NA_letter
#' @template three-dots
#'
#' @return An \code{\link[=sq-class]{sq}} object of the same type as input
#' object but each sequence is reversed.
#' 
#' @details
#' \code{reverse()} function reverses each sequence in supplied \code{sq} object
#' (e.q. transforms \code{"MIAANYTWIL"} to \code{"LIWTYNAAIM"}). This operation
#' does not alter the type of the input object nor its alphabet.
#'
#' @examples
#' # Creating objects to work on:
#' sq_ami <- sq(c("MIAANYTWIL","TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"),
#'              alphabet = "ami_ext")
#' sq_dna <- sq(c("ATGCAGGA", "GACCGAACGAN", ""), alphabet = "dna_ext")
#' sq_unt <- sq(c("ATGCAGGA?", "TGACGAGCTTA", "", "TIAALGNIIYRAIE"))
#'
#' # Reversing sequences:
#' reverse(sq_ami)
#' reverse(sq_dna)
#' reverse(sq_unt)
#'
#' @family order_functions
#' @export
reverse <- function(x, ...)
  UseMethod("reverse")

#' @export
reverse.default <- function(x, ...)
  stop("'reverse' isn't implemented for this type of object; maybe you wanted to use 'rev'?", call. = FALSE)

#' @rdname reverse
#' @export
reverse.sq <- function(x, ...,
                       NA_letter = getOption("tidysq_NA_letter")) {
  assert_string(NA_letter, min.chars = 1)
  
  CPP_reverse(x, NA_letter)
}
