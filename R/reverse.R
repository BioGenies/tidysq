#' Reverse sequence
#' 
#' @description Reverse given list of sequences.
#' 
#' @template x
#' @template NA_letter
#' @template three-dots
#' 
#' @details The \code{reverse} function reverses each sequence in supplied 
#' \code{\link[=sq-class]{sq}} object (e.q. transforms "MIAANYTWIL" to "LIWTYNAAIM").
#' Empty sequences are left with no effect. This operation does not change 
#' the type of the input object nor its alphabet.
#' 
#' Since the function \code{reverse} returns a \code{\link[=sq-class]{sq}} object, the
#' \code{\link[=sq-print]{print}} function is implicitly called.
#' 
#' @return A \code{\link[=sq-class]{sq}} object of the same type as input object but
#' each sequence is reversed.
#' 
#' @seealso \code{\link[=sq-class]{sq}} \code{\link{clean}} \code{\link{sq-print}}
#' 
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
