#' Convert DNA or RNA into proteins using genetic code
#' 
#' @description This function allows the user to input DNA or RNA sequences and
#' acquire sequences of corresponding proteins, where correspondence is encoded
#' in specified table.
#'
#' @param x [\code{sq_dna_bsc} || \code{sq_rna_bsc}]\cr
#'  An object this function is applied to.
#' @param table [\code{integer(1)}]\cr
#'  The number of translation table used, as specified
#' \href{https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi}{here}.
#' @template NA_letter
#' @template three-dots
#' 
#' @return An object of \code{\link[=sq-class]{class sq}} with \strong{ami_bsc}
#' type.
#' 
#' @details
#' DNA and RNA sequences use combinations of three consecutive nucleic acids to
#' encode one of 22 amino acids. This encoding is called "genetic code".
#' 
#' \code{translate()} first splits passed DNA or RNA sequences into
#' three-letter chunks. Then searches the codon table for the entry where the
#' key is equal to the current chunk and the value is one letter that encodes
#' the corresponding protein. These resulting letters are then pasted into one
#' sequence for each input sequence.
#' 
#' Due to how the tables works, \code{translate()} does not support inputting
#' sequences with extended alphabets, as ambiguous letters in most cases cannot
#' be translated into exactly one protein.
#' 
#' Moreover, this function raises an error whenever input sequence contain
#' either "\code{-}" or \code{NA} value.
#'
#' @examples
#' sq_dna <- sq(c("TACTGGGCATGA", "CAGGTC", "TAGTCCTAG"), alphabet = "dna_bsc")
#' translate(sq_dna)
#'
#' @family bio_functions
#' @seealso \code{\link{remove_ambiguous}}, \code{\link{substitute_letters}} and
#' \code{\link{typify}} for necessary actions before using \code{translate()}
#' @export
translate <- function(x, table = 1, ...) {
  assert_int(table)
  assert_choice(table, c(1:16, 21:26, 29, 30, 33))
  
  UseMethod("translate")
}

#' @export
translate.default <- function(x, table = 1, ...)
  stop("cannot translate something that is neither basic DNA nor RNA sequence", call. = FALSE)

#' @rdname translate
#' @export
translate.sq_dna_bsc <- function(x, table = 1, ...,
                                 NA_letter = getOption("tidysq_NA_letter")) {
  assert_string(NA_letter, min.chars = 1)
  
  CPP_translate(x, table, NA_letter)
}

#' @rdname translate
#' @export
translate.sq_rna_bsc <- translate.sq_dna_bsc
