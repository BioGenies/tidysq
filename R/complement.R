#' Create complement sequence from dnasq or rnasq object 
#' 
#' @description Creates the complementary sequence from a given RNA or DNA 
#'  nucleotide sequence. The function differentiates between RNA and DNA sequences. 
#' 
#' @param x a \code{\link[=sq-class]{sq}} object of type \strong{dna} or \strong{rna}.
#'
#' @return \code{sq} object of the same type as input \code{dnasq} (\strong{dna})
#' or \code{rnasq} (\strong{rna}) but built of complementary nucleotides to entered
#' sequence.
#' 
#' @details This function allows to get complement sequence, which is created by 
#' matching elements (nucleotides) with complementary to input dnasq or rnasq object.
#' Whether 'U' (uracil) or 'T' (thymine) is used depends on the class of the sq object.
#' 
#' Functions \code{complement_dna} and \code{complement_rna} are provided as a safe
#' way of limiting classes \code{complement} function is used on.
#' 
#' @seealso \code{\link[=sq-class]{sq}}
#' @export
complement <- function(x, ...)
  UseMethod("complement")

#' @export
complement.default <- function(x, ...)
  stop("method 'complement' isn't implemented for this type of object", call. = FALSE)

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
