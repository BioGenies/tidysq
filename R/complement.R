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
  
  inds_fun <- match(dict[alph], alph) - 1
  ret <- pack(lapply(ret, function(s) inds_fun[s + 1]), alph)
  
  vec_restore(ret, x)
}

#' @export
complement.sq_rna_bsc <- function(x) {
  alph <- alphabet(x)
  ret <- unpack(x, "INTS")
  
  dict <- c(G = "C", C = "G", U = "A", A = "U", `-` = "-")
  
  inds_fun <- match(dict[alph], alph) - 1
  ret <- pack(lapply(ret, function(s) inds_fun[s + 1]), alph)
  
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
