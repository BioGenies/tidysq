#' @export
remove_ambiguous <- function(x, by_letter = FALSE, ...) {
  assert_flag(by_letter)
  UseMethod("remove_ambiguous")
}

#' @export
remove_ambiguous.default <- function(x, by_letter = FALSE, ...)
  stop("ambiguous letters are not defined in the context of this class", call. = FALSE)

#' @export
remove_ambiguous.sq_dna_bsc <- function(x, by_letter = FALSE, ...) x

#' @export
remove_ambiguous.sq_dna_ext <- function(x, by_letter = FALSE, ...)
  CPP_remove_ambiguous(x, by_letter, get_standard_alphabet("sq_dna_bsc"))

#' @export
remove_ambiguous.sq_rna_bsc <- function(x, by_letter = FALSE, ...) x

#' @export
remove_ambiguous.sq_rna_ext <- function(x, by_letter = FALSE, ...)
  CPP_remove_ambiguous(x, by_letter, get_standard_alphabet("sq_rna_bsc"))

#' @export
remove_ambiguous.sq_ami_bsc <- function(x, by_letter = FALSE, ...) x

#' @export
remove_ambiguous.sq_ami_ext <- function(x, by_letter = FALSE, ...)
  CPP_remove_ambiguous(x, by_letter, get_standard_alphabet("sq_ami_bsc"))
