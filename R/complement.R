#' Create complement sequence from dnasq or rnasq object 
#' 
#' @description Creates the complementary sequence from a given RNA or DNA 
#'  nucleotide sequence. The function differentiates between RNA and DNA sequences. 
#' 
#' @param x a \code{\link{sq}} object of type \strong{dna} or \strong{rna}.
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
#' @examples 
#' # Creating objects dna_sequence (with DNA sequences) and rna_sequence 
#' # (with RNA sequences) to work on:
#' 
#' dna_sequence <- construct_sq(c("ACTGCTG", "CTTAGA", "CCCT", "CTGAATGT"),
#'                              type = "dna")
#' rna_sequence <- construct_sq(c("ACUGCUG", "CUUAGA", "CCCU", "CUGAAUGU"),
#'                              type = "rna")
#' 
#' 
#' # In the following the complement function is used to make a PCR (Polymerase Chain Reaction)
#' # primers. Every sequence will be rewritten to its complementary equivalent as 
#' # following example: AAATTTGGG to TTTAAACCC.
#'  
#' # creating complementary sequences with the basic function and using wrappers:
#' complement(dna_sequence)
#' complement_dna(dna_sequence)
#' 
#' complement(rna_sequence)
#' complement_rna(rna_sequence)
#' 
#' # Each sequence from dna_sequence and rna_sequence object have now an own
#' # complementary equivalent, which can be helpful during constructing PCR primers.
#' 
#' @seealso \code{\link{sq}}
#' @export
complement <- function(x)
  UseMethod("complement")

#' @export
complement.default <- function(x)
  stop("method 'complement' isn't implemented for this type of object", call. = FALSE)

#' @export
complement.sq_dna_bsc <- function(x) {
  alph <- alphabet(x)
  ret <- unpack(x, "INTS")
  
  dict <- c(G = "C", C = "G", T = "A", A = "T", `-` = "-")
  
  inds_fun <- match(dict[alph], alph)
  names(inds_fun) <- as.character(1:length(alph))
  ret <- lapply(ret, function(s) structure(pack(inds_fun[s], "INTS"),
                                           original_length = attr(s, "original_length")))
  
  vec_restore(ret, x)
}

#' @export
complement.sq_rna_bsc <- function(x) {
  alph <- alphabet(x)
  ret <- unpack(x, "INTS")
  
  dict <- c(G = "C", C = "G", U = "A", A = "U", `-` = "-")
  
  inds_fun <- match(dict[alph], alph)
  names(inds_fun) <- as.character(1:length(alph))
  ret <- lapply(ret, function(s) structure(pack(inds_fun[s], "INTS"),
                                           original_length = attr(s, "original_length")))
  
  vec_restore(ret, x)
}

#' @rdname complement
#' @export
complement_dna <- function(x)
  UseMethod("complement_dna")

#' @export
complement_dna.default <- function(x)
  stop("method 'complement_dna' isn't implemented for this type of object", call. = FALSE)

#' @export
complement_dna.sq_dna_bsc <- complement.sq_dna_bsc

#' @rdname complement
#' @export
complement_rna <- function(x)
  UseMethod("complement_rna")

#' @export
complement_rna.default <- function(x)
  stop("method 'complement_rna' isn't implemented for this type of object", call. = FALSE)

#' @export
complement_rna.sq_rna_bsc <- complement.sq_rna_bsc
