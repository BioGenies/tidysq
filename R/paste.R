#' @export
paste <- function(...)
  UseMethod("paste")

#' @export
paste.default <- function(...) {
  base::paste(...)
}

#' Paste sequences in string-like fashion
#'
#' @description Joins multiple vectors of sequences into one vector.
#'
#' @param ... [\code{sq}]\cr
#'  Sequences to paste together.
#' @template NA_letter
#'
#' @return \code{\link[=sq-class]{sq}} object of common type of input objects.
#' Common type is determined in the same process as for
#' \code{\link[=sq-concatenate]{c.sq}()}.
#'
#' @details
#' \code{paste()} joins sequences in the same way as it does with strings.
#' All \code{sq} objects must have the same length, that is, contain the same
#' number of sequences. An exception is made for scalar (length 1) \code{sq}
#' objects, which are replicated instead.
#'
#' @examples
#' # Creating objects to work on:
#' sq_dna_1 <- sq(c("TTCAGGGCTAG", "CGATTGC", "CAGTTTA"),
#'                alphabet = "dna_bsc")
#' sq_dna_2 <- sq(c("ATCTTGAAG", "CATATGCGCTA", "ACGTGTCGA"),
#'                alphabet = "dna_bsc")
#' sq_unt_1 <- sq(c("ATGCAGGA?", "TGACGAGCTTA", "", "TIAALGNIIYRAIE"))
#' sq_unt_2 <- sq(c("OVNU!!OK!!J", "GOK!MI!N!BB!", "DPOFIN!!", "??!?"))
#'
#' # Pasting sequences:
#' collapse(sq_dna_1, sq_dna_2)
#' collapse(sq_unt_1, sq_unt_2)
#' collapse(sq_dna_2, sq_unt_2, sq_dna_1)
#'
#' @family order_functions
#' @name paste
#' @export
paste.sq <- function(...,
                     NA_letter = getOption("tidysq_NA_letter")) {
  # Throws error when there is no common size
  vec_size_common(...)
  assert_string(NA_letter, min.chars = 1)

  CPP_paste(vec_cast_common(...), NA_letter)
}
