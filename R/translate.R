#' Convert DNA or RNA into proteins using genetic code
#' 
#' @description This function allows the user to input DNA or RNA sequences and
#' acquire sequences of corresponding proteins, where correspondence is encoded
#' in specified table.
#' 
#' @param x an object of class \code{\link{sq}} with either \strong{dna} or
#' \strong{rna} type
#' @param table integer number of translation table used, as specified
#' \href{https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi}{here}
#' 
#' @return An object of \code{\link[=sq]{class sq}} with \strong{ami} type and
#' \strong{cln} subtype.
#' 
#' @details
#' DNA and RNA sequences use combinations of three consecutive nucleic acids to
#' encode one of 22 amino acids. This encoding is called "genetic code".
#' 
#' \code{translate()} first splits passed DNA or RNA sequences into
#' three-letter chunks. Then searches the codon table for the entry where the
#' key is equal to the current chunk and the value is one letter that encodes
#' the corresponding protein. These resulting letters are then pasted into one
#' sequence for each input sequence and objectifized into an \code{sq} object.
#' 
#' Due to how the table works, \code{translate()} does not support inputting
#' uncleaned sequences, as ambiguous letters mostly cannot be translated into
#' exactly one protein.
#' 
#' Moreover, behaviour of this function is undefined whenever cleaned sequence
#' contain either \code{-} or \code{NA} sign.
#' 
#' @seealso \code{\link{clean}}, \code{\link{substitute_letters}} and
#' \code{\link{typify}} for necessary actions before using \code{translate()}
#' 
#' @examples
#' sq_dna <- construct_sq_dna(c("TACTGGGCATGA", "CAGGTC", "TAGTCCTAG"),
#'                            is_clean = TRUE)
#' translate(sq_dna)
#' 
#' @export
translate <- function(x, table = 1,
                      NA_letter = getOption("tidysq_NA_letter"), ...) {
  assert_int(table)
  assert_choice(table, c(1:16, 21:31, 33))
  assert_string(NA_letter)
  
  UseMethod("translate")
}

#' @export
translate.default <- function(x, table = 1,
                              NA_letter = getOption("tidysq_NA_letter"), ...)
  stop("cannot translate something that is neither basic DNA nor RNA sequence", call. = FALSE)

#' @export
translate.sq_dna_bsc <- function(x, table = 1,
                                 NA_letter = getOption("tidysq_NA_letter"), ...,
                                 interpret_as_stop = FALSE) {
  assert_flag(interpret_as_stop)
  CPP_translate(x, table, NA_letter, interpret_as_stop)
}

#' @export
translate.sq_rna_bsc <- translate.sq_dna_bsc
