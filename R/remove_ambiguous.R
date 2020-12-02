#' Remove sequences that contain ambiguous elements
#'
#' @param x [\code{sq_dna_bsc} || \code{sq_rna_bsc} || \code{sq_dna_ext} ||
#' \code{sq_rna_ext} || \code{sq_ami_bsc} || \code{sq_ami_ext}]\cr
#'  An object this function is applied to.
#' @template by_letter
#' @template NA_letter
#' @template three-dots
#'
#' @family cleaning_functions
#' @export
remove_ambiguous <- function(x, by_letter = FALSE, ...) {
  assert_flag(by_letter)
  UseMethod("remove_ambiguous")
}

#' @export
remove_ambiguous.default <- function(x, by_letter = FALSE, ...)
  stop("ambiguous letters are not defined in the context of this class", call. = FALSE)

#' @export
remove_ambiguous.sq_dna_bsc <- function(x, by_letter = FALSE, ...,
                                        NA_letter = getOption("tidysq_NA_letter")) x

#' @export
remove_ambiguous.sq_dna_ext <- function(x, by_letter = FALSE, ...,
                                        NA_letter = getOption("tidysq_NA_letter")) {
  assert_string(NA_letter, min.chars = 1)
  
  CPP_remove_ambiguous(x, by_letter, NA_letter)
}

#' @export
remove_ambiguous.sq_rna_bsc <- function(x, by_letter = FALSE, ...,
                                        NA_letter = getOption("tidysq_NA_letter")) x

#' @export
remove_ambiguous.sq_rna_ext <- function(x, by_letter = FALSE, ...,
                                        NA_letter = getOption("tidysq_NA_letter")) {
  assert_string(NA_letter, min.chars = 1)
  
  CPP_remove_ambiguous(x, by_letter, NA_letter)
}

#' @export
remove_ambiguous.sq_ami_bsc <- function(x, by_letter = FALSE, ...,
                                        NA_letter = getOption("tidysq_NA_letter")) x

#' @export
remove_ambiguous.sq_ami_ext <- function(x, by_letter = FALSE, ...,
                                        NA_letter = getOption("tidysq_NA_letter")) {
  assert_string(NA_letter, min.chars = 1)
  
  CPP_remove_ambiguous(x, by_letter, NA_letter)
}
