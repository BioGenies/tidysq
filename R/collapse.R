#' Collapse multiple sequences into one
#'
#' @description Joins sequences from a vector into a single sequence. Sequence
#' type remains unchanged.
#'
#' @template x
#' @template NA_letter
#' @template three-dots
#'
#' @return \code{\link[=sq-class]{sq}} object of the same type as input but with
#' exactly one sequence.
#'
#' @details
#' \code{collapse()} joins sequences from supplied \code{sq} object in the same
#' order as they appear in said vector. That is, if there are three sequences
#' AGGCT, ATCCGT and GAACGT, then resulting sequence will be AGGCTATCCGTGAACGT.
#' This operation does not alter the type of the input object nor its alphabet.
#'
#' @examples
#' # Creating objects to work on:
#' sq_ami <- sq(c("MIAANYTWIL","TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"),
#'              alphabet = "ami_ext")
#' sq_dna <- sq(c("ATGCAGGA", "GACCGAACGAN", ""), alphabet = "dna_ext")
#' sq_unt <- sq(c("ATGCAGGA?", "TGACGAGCTTA", "", "TIAALGNIIYRAIE"))
#'
#' # Collapsing sequences:
#' collapse(sq_ami)
#' collapse(sq_dna)
#' collapse(sq_unt)
#'
#' # Empty sq objects are collapsed as well (into empty string - ""):
#' sq_empty <- sq(character(), alphabet = "rna_bsc")
#' collapse(sq_empty)
#'
#' @family order_functions
#' @export
collapse <- function(x, ...)
  UseMethod("collapse")

#' @export
collapse.default <- function(x, ...)
  stop_no_method(collapse, x)

#' @rdname collapse
#' @export
collapse.sq <- function(x, ...,
                        NA_letter = getOption("tidysq_NA_letter")) {
  assert_string(NA_letter, min.chars = 1)

  CPP_collapse(x, NA_letter)
}
