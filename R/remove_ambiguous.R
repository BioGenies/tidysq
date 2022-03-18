#' Remove sequences that contain ambiguous elements
#'
#' @description This function replaces sequences with ambiguous elements by
#' empty (\code{NULL}) sequences or removes ambiguous elements from sequences
#' in an \code{sq} object.
#'
#' @param x [\code{sq_dna_bsc || sq_rna_bsc || sq_dna_ext || sq_rna_ext || sq_ami_bsc || sq_ami_ext}]\cr
#'  An object this function is applied to.
#' @template by_letter
#' @template NA_letter
#' @template three-dots
#'
#' @return An \code{\link[=sq-class]{sq}} object with the \strong{_bsc}
#' version of inputted type.
#'
#' @details
#' Biological sequences, whether of DNA, RNA or amino acid elements, are not
#' always exactly determined. Sometimes the only information the user has about
#' an element is that it's one of given set of possible elements. In this case
#' the element is described with one of special letters, here called
#' \strong{ambiguous}.
#'
#' The inclusion of these letters is the difference between extended and basic
#' alphabets (and, conversely, types). For amino acid alphabet these letters
#' are: B, J, O, U, X, Z; whereas for DNA and RNA: W, S, M, K, R, Y, B, D, H, V,
#' N.
#'
#' \code{remove_ambiguous()} is used to create sequences without any of the
#' elements above. Depending on value of \code{by_letter} argument, the function
#' either replaces "ambiguous" sequences with empty sequences (if
#' \code{by_letter} is equal to \code{TRUE}) or shortens original sequence by
#' retaining only unambiguous letters (if opposite is true).
#'
#' @examples
#' # Creating objects to work on:
#' sq_ami <- sq(c("MIAANYTWIL","TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"),
#'              alphabet = "ami_ext")
#' sq_dna <- sq(c("ATGCAGGA", "GACCGAACGAN", "TGACGAGCTTA", "ACTNNAGCN"),
#'              alphabet = "dna_ext")
#'
#' # Removing whole sequences with ambiguous elements:
#' remove_ambiguous(sq_ami)
#' remove_ambiguous(sq_dna)
#'
#' # Removing ambiguous elements from sequences:
#' remove_ambiguous(sq_ami, by_letter = TRUE)
#' remove_ambiguous(sq_dna, by_letter = TRUE)
#'
#' # Analysis of the result
#' sq_clean <- remove_ambiguous(sq_ami)
#' is_empty_sq(sq_clean)
#' sq_type(sq_clean)
#'
#' @family cleaning_functions
#' @export
remove_ambiguous <- function(x, by_letter = FALSE, ...) {
  assert_flag(by_letter)
  
  UseMethod("remove_ambiguous")
}

#' @export
remove_ambiguous.default <- function(x, by_letter = FALSE, ...)
  stop_no_method(remove_ambiguous, x)

#' @rdname remove_ambiguous
#' @export
remove_ambiguous.sq <- function(x, by_letter = FALSE, ...,
                                NA_letter = getOption("tidysq_NA_letter")) {
  assert_string(NA_letter, min.chars = 1)

  CPP_remove_ambiguous(x, by_letter, NA_letter)
}
