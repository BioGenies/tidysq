#' Subset sequences from sq objects
#' 
#' @description Extracts a defined range of elements from all sequences.
#' 
#' @template x
#' @param indices [\code{integer}]\cr
#'  Indices to extract from each sequence. The function follows the normal R
#'  conventions for indexing vectors, including negative indices.
#' @template three-dots
#' @template NA_letter
#' @template on_warning
#' 
#' @return \code{\link[=sq-class]{sq}} object of the same type as input
#' \code{sq}, where each element is a subsequence created by indexing
#' corresponding sequence from input \code{sq} object with input indices.
#' 
#' @details
#' \code{bite} function allows user to access specific elements from multiple
#' sequences at once.
#'
#' By passing positive indices the user can choose, which elements they want
#' from each sequence. If a sequence is shorter than an index, then \code{NA}
#' value is inserted into the result in this place and a warning is issued.
#' The user can specify behavior of R in this case by specifying
#' \code{on_warning} parameter.
#'
#' Negative indices are supported as well. Their interpretation is "to select
#' all elements except those on positions specified by these negative indices".
#' This means that e.g. \code{c(-1, -3, -5)} vector will be used to bite all
#' sequence elements except the first, the third and the fifth. If a sequence
#' is shorter than any index, then nothing happens, as it's physically
#' impossible to extract an element at said index.
#'
#' As per normal R convention, it isn't accepted to mix positive and negative
#' indices, because there is no good interpretation possible for that.
#'
#' @examples
#' # Creating objects to work on:
#' sq_dna <- sq(c("ATGCAGGA", "GACCGNBAACGAN", "TGACGAGCTTA"),
#'              alphabet = "dna_bsc")
#' sq_ami <- sq(c("MIAANYTWIL","TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"),
#'              alphabet = "ami_ext")
#' sq_unt <- sq(c("ATGCAGGA?", "TGACGAGCTTA", "", "TIAALGNIIYRAIE"))
#'
#' # Extracting first five letters:
#' bite(sq_dna, 1:5)
#'
#' # If a sequence is shorter than 5, then NA is introduced:
#' bite(sq_unt, 1:5)
#'
#' # Selecting fourth, seventh and fourth again letter:
#' bite(sq_ami, c(4, 7, 4))
#'
#' # Selecting all letters except first four:
#' bite(sq_dna, -1:-4)
#'
#' @family order_functions
#' @seealso \code{\link{remove_na}}
#' @export
bite <- function(x, indices, ...)
  UseMethod("bite")

#' @export
bite.default <- function(x, indices, ...)
  stop("method 'bite()' isn't implemented for this type of object", call. = FALSE)

#' @rdname bite
#' @export
bite.sq <- function(x, indices, ...,
                    NA_letter = getOption("tidysq_NA_letter"),
                    on_warning = getOption("tidysq_on_warning")) {
  assert_string(NA_letter, min.chars = 1)
  assert_warning_handling(on_warning)
  assert_integerish(indices, any.missing = FALSE, null.ok = TRUE)
  
  CPP_bite(x, indices, NA_letter, on_warning)
}
