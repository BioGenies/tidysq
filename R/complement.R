#' Create complement sequence from dnasq or rnasq object 
#' 
#' @description Creates the complementary sequence from a given RNA or DNA 
#' sequence. The function keeps the type of sequence intact.
#' 
#' @param x [\code{sq_dna_bsc} || \code{sq_rna_bsc} || \code{sq_dna_ext} ||
#' \code{sq_rna_ext}]\cr
#'  An object this function is applied to.
#' @template NA_letter
#' @template three-dots
#'
#' @return \code{\link[=sq-class]{sq}} object of the same type as input but
#' built of nucleotides complementary to those in the entered sequences.
#' 
#' @details
#' This function matches elements of sequence to their complementary letters.
#' For unambiguous letters, "\code{C}" is matched with "\code{G}" and "\code{A}"
#' is matched with either "\code{T}" (thymine) or "\code{U}" (uracil), depending
#' on whether input is of \strong{dna} or \strong{rna} type.
#'
#' Ambiguous letters are matched as well, for example "\code{N}" (any
#' nucleotide) is matched with itself, while "\code{B}" (not alanine) is matched
#' with "\code{V}" (not thymine/uracil).
#'
#' @examples
#' # Creating DNA and RNA sequences to work on:
#' sq_dna <- sq(c("ACTGCTG", "CTTAGA", "CCCT", "CTGAATGT"),
#'              alphabet = "dna_bsc")
#' sq_rna <- sq(c("BRAUDUG", "URKKBKUCA", "ANKRUGBNNG", "YYAUNAAAG"),
#'              alphabet = "rna_ext")
#'
#' # Here complement() function is used to make PCR (Polymerase Chain Reaction)
#' # primers. Every sequence is rewritten to its complementary equivalent as
#' # in the following example: AAATTTGGG -> TTTAAACCC.
#'
#' complement(sq_dna)
#' complement(sq_rna)
#'
#' # Each sequence have now a complementary equivalent, which can be helpful
#' # during constructing PCR primers.
#'
#' @family bio_functions
#' @seealso \code{\link[=sq-class]{sq}}
#' @export
complement <- function(x, ...)
  UseMethod("complement")

#' @export
complement.default <- function(x, ...)
  stop_no_method(complement, x)

#' @rdname complement
#' @export
complement.sq_dna_bsc <- function(x, ...,
                                  NA_letter = getOption("tidysq_NA_letter")) {
  assert_string(NA_letter, min.chars = 1)
  
  CPP_complement(x, NA_letter)
}

#' @rdname complement
#' @export
complement.sq_dna_ext <- complement.sq_dna_bsc

#' @rdname complement
#' @export
complement.sq_rna_bsc <- complement.sq_dna_bsc

#' @rdname complement
#' @export
complement.sq_rna_ext <- complement.sq_dna_bsc
