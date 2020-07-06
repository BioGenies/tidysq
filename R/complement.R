#' Create complement sequence from nucsq object 
#' 
#' @description Creates the complementary sequence from a given RNA or DNA 
#'  nucleotide sequence. The function differentiates between RNA and DNA sequences. 
#' 
#' @param nucsq a \code{\link{sq}} object of type \strong{nuc}.
#' @param is_dna a \code{\link{logical}} indicating if entered sequence is DNA or RNA. If 
#' \code{is_dna} is \code{TRUE}, sequences cannot contain 'U' (uracil);
#' if is \code{FALSE}, sequences cannot contain 'T' (thymine). If 
#' \code{NULL} (default) the sequence type is not specified and is guessed (see details below).
#'
#' @return \code{sq} object of the same type as input \code{nucsq} (\strong{nuc}) but 
#' built of complementary nucleotides to entered sequence.
#' 
#' @details This function allows to get complement sequence, which is created by 
#' matching elements (nucleotides) with complementary to input nucsq object. If 
#' \code{is_dna = TRUE} entered sequence is DNA. If \code{is_dna = FALSE} entered 
#' sequence is RNA. By default, the sequence type is not specified, and the function
#' tries to guess, which type of sequence was entered. Sequences containing 'U' without 
#' 'T' will be set to the type RNA. If a sequence contains 'T' (thymine) without 
#' 'U' (uracil), the type is set to DNA. An error is displayed if both 'T' and 'U' 
#' are present in the sequence or if the sequence contains only 'A' (adenine). If the 
#' sequence does not contain 'T' or 'U' or the logical specification is wrong 
#' (i.e., if the sequence contains 'U' and the logical specification is set to DNA), 
#' an error will also be returned.
#' 
#' Functions \code{complement_dna} and \code{complement_rna} are wrappers around 
#' \code{complement} with specified \code{is_dna} parameter - accordingly TRUE or FALSE. 
#' 
#' @examples 
#' # Creating objects nuc_dna_sequence (with DNA sequences) and nuc_rna_sequence 
#' # (with RNA sequences) to work on:
#' 
#' nuc_dna_sequence <- construct_sq(c("ACTGCTG", "CTTAGA", 
#'                                    "CCCT", "CTGAATGT"), type = "nuc")
#' nuc_rna_sequence <- construct_sq(c("ACUGCUG", "CUUAGA", 
#'                                    "CCCU", "CUGAAUGU"), type = "nuc")
#'                                    
#' 
#' # In the following the complement function is used to make a PCR (Polymerase Chain Reaction)
#' # primers. Every sequence will be rewritten to its complementary equivalent as 
#' # following example: AAATTTGGG to TTTAAACCC.
#'  
#' # creating complementary sequences with a specified sequence type (and using wrappers):
#' complement(nuc_dna_sequence, is_dna = TRUE)
#' complement_dna(nuc_dna_sequence)
#' 
#' complement(nuc_rna_sequence, is_dna = FALSE)
#' complement_rna(nuc_rna_sequence)
#' 
#' # creating complementary sequences without defined sequence type:
#' complement(nuc_dna_sequence)
#' complement(nuc_rna_sequence)
#' 
#' # Each sequence from nuc_dna_sequence and nuc_rna_sequence object have now an  
#' # own complementary equivalent, which can be helpful during constructing PCR primers.
#' 
#' @seealso \code{\link{sq}}
#' @export
complement <- function(sq) {
  UseMethod("complement")
}

#' @export
complement.default <- function(sq) {
  stop("method 'complement' isn't implemented for this type of object")
}

#' @export
complement.dnasq <- function(sq) {
  .validate_sq(sq, "dna")
  
  .check_is_clean(sq, "'dnasq'")
  alph <- .get_alph(sq)
  alph_size <- .get_alph_size(alph)
  sq <- .unpack_from_sq(sq, "int")
  
  dict <- c(G = "C", C = "G", T = "A", A = "T", `-` = "-")
  
  inds_fun <- match(dict[alph], alph)
  names(inds_fun) <- as.character(1:length(alph))
  ret <- lapply(sq, function(s)  C_pack_ints(inds_fun[s], alph_size))
  
  ret <- .set_alph(ret, alph)
  .set_class(ret, "dna", TRUE)
}

#' @export
complement.rnasq <- function(sq) {
  .validate_sq(sq, "rna")
  
  .check_is_clean(sq, "'rnasq'")
  alph <- .get_alph(sq)
  alph_size <- .get_alph_size(alph)
  sq <- .unpack_from_sq(sq, "int")
  
  dict <- c(G = "C", C = "G", U = "A", A = "U", `-` = "-")
  
  inds_fun <- match(dict[alph], alph)
  names(inds_fun) <- as.character(1:length(alph))
  ret <- lapply(sq, function(s)  C_pack_ints(inds_fun[s], alph_size))
  
  ret <- .set_alph(ret, alph)
  .set_class(ret, "rna", TRUE)
}

#' @rdname complement
#' @export
complement_dna <- function(sq) {
  UseMethod("complement")
}

#' @export
complement_dna.default <- function(sq) {
  stop("method 'complement_dna' isn't implemented for this type of object")
}

#' @export
complement_dna.dnasq <- function(sq) {
  complement.dnasq(sq)
}

#' @rdname complement
#' @export
complement_rna <- function(sq) {
  UseMethod("complement")
}

#' @export
complement_rna.default <- function(sq) {
  stop("method 'complement_rna' isn't implemented for this type of object")
}

#' @export
complement_rna.dnasq <- function(sq) {
  complement.rnasq(sq)
}
