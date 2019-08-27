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
#' @details This functions allow obtaining complement sequence which is created by 
#' matching elements (nucleotides) with complementarity to input \code{nucsq} object. In case
#' of \code{complement}: If \code{is_dna = TRUE} entered sequence is DNA. 
#' If \code{is_dna = FALSE} entered 
#' sequence is RNA. By default the sequence type is not specified and the function
#' tries to guess which type of sequence was entered. If sequence contain 'U' without 
#' 'T' the type will be set to RNA. If a sequence contains 'T' (thymine) without 
#' 'U' (uracil), the type is set to DNA. An error is displayed if both 'T' and 'U' 
#' are present in the sequence or if the sequence contains only 'A' (adenine). If the 
#' sequence does not contain 'T' or 'U' or the logical specification is wrong 
#' (i.e. if the sequence contains 'U' and the logical specification is set to DNA), 
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
#' # creating complementary sequences with defined sequence type (and using wrappers):
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
#' # Each sequence from  nuc_dna_sequence and nuc_rna_sequence object have now own 
#' # complementary equivalent, which can be helpful during constructing PCR primers.
#'  
#'   
#' 
#' @seealso \code{\link{sq}}
#' 
#' @export
complement <- function(nucsq, is_dna = NULL) {
  validate_sq(nucsq, "nuc")
  
  .check_is_clean(nucsq, "'nucsq'")
  alph <- .get_alph(nucsq)
  alph_size <- .get_alph_size(alph)
  sq <- .debitify_sq(nucsq, "int")
  
  has_U <- any(unlist(sq) == match("U", alph))
  has_T <- any(unlist(sq) == match("T", alph))
  has_A <- any(unlist(sq) == match("A", alph))
  .check_has_both_UT(has_U, has_T)
  dict <- c(G = "C", C = "G", T = "A", U = "A", A = "T", `-` = "-")
  
  if (is.null(is_dna)) {
    .check_has_A_no_UT(has_U, has_T, has_A) 
    if (has_U) dict["A"] <- "U"
  } else {
    .check_logical(is_dna, "'is_dna'", allow_null = TRUE)
    .check_is_dna_matches(is_dna, has_U, has_T) 
    if (!is_dna) dict["A"] <- "U"
  }
  
  inds_fun <- match(dict[alph], alph)
  names(inds_fun) <- as.character(1:length(alph))
  ret <- lapply(sq, function(s) pack_ints(inds_fun[s], alph_size))
  
  ret <- .set_alph(ret, alph)
  .set_class(ret, "nuc", TRUE)
}

#' @rdname complement
#' @export
complement_dna <- function(nucsq) {
  complement(nucsq, is_dna = TRUE)
}

#' @rdname complement
#' @export
complement_rna <- function(nucsq) {
  complement(nucsq, is_dna = FALSE)
}