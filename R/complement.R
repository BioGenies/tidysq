#' Create complement sequence from nucsq object 
#' 
#' @description Creates the complementary sequence from a given RNA or DNA 
#'  nucleotide sequence. The function differentiates between RNA and DNA sequences. 
#' 
#' @param nucsq \code{\link{nucsq}} object of type 'nuc'.
#' @param is_dna logical indicating if entered sequence is DNA or RNA. If 
#' 'is_dna' is \code{TRUE}, sequences cannot contain 'U' (uracil);
#' if is \code{FALSE}, sequences cannot contain 'T' (thymine). If 
#' \code{NULL} (default) the sequence type is not specified.
#'
#' @return \code{\link{nucsq}} object of the same type as input nucsq but 
#' built of complementary nucleotides to entered sequence.
#' 
#' @details This function allows to get complement sequence, which is created by 
#' matching elements (nucleotides) with complementary to input nucsq object. If 
#' \code{is_dna = TRUE} entered sequence is DNA. If \code{is_dna = FALSE} entered 
#' sequence is RNA. By default, the sequence type is not specified, and the function
#' tries to guess, which type of sequence was entered. Sequences containing 'U' without 
#' 'T' will be set to the tpye RNA. If a sequence contains 'T' (thymine) without 
#' 'U' (uracil), the type is set to dna. An error is displayed if both 'T' and 'U' 
#' are present in the sequence or if the sequence contains only 'A' (adenine). If the 
#' sequence does not contain 'T' or 'U' or the logical specification is wrong 
#' (i.e., if the sequence contains 'U' and the logical specification is set to DNA), 
#' an error will also be returned.
#' Both RNA and DNA sequences can be rewritten to complementary sequence. 
#' 
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
#' # Sequences are now ready to be used as matrices for creating complementary sequences.
#' # Get an overview of the sequences:
#' summary(nuc_dna_sequence)
#' summary(nuc_rna_sequence)
#' 
#' # In the following the complement function is used to make a PCR (Polymerase Chain Reaction)
#' # primers. Every sequence will be rewritten to its complementary equivalent as 
#' # following example: AAATTTGGG to TTTAAACCC.
#'  
#' # creating complementary sequences with defined sequence type:
#' complement(nuc_dna_sequence, is_dna = TRUE)
#' complement(nuc_rna_sequence, is_dna = FALSE)
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
#' @seealso sq clnsq nucsq
#' 
#' @export
complement <- function(nucsq, is_dna = NULL) {
  validate_sq(nucsq, "nuc")
  if (!.is_cleaned(nucsq)) {
    stop("'nucsq' needs to be cleaned")
  }
  alph <- .get_alph(nucsq)
  alph_size <- .get_alph_size(alph)
  sq <- .debitify_sq(nucsq, "int")
  
  has_U <- any(unlist(sq) == match("U", alph))
  has_T <- any(unlist(sq) == match("T", alph))
  has_A <- any(unlist(sq) == match("A", alph))
  
  if (has_U && has_T) {
    stop("'nucsq' sequences contains both 'U' and 'T' letters - should contain only one of them")
  }
  
  dict <- c(G = "C", C = "G", T = "A", U = "A", A = "T", `-` = "-")
  if (is.null(is_dna)) {
    if (!has_U && !has_T && has_A) {
      stop("'nucsq' sequences contains 'A' elements, but does not contain 'T' nor 'U' - unable to guess if it's dna or rna")
    } 
    if (has_U) dict["A"] <- "U"
  } else {
    if (!is.logical(is_dna) ||
        is.na(is_dna) ||
        (length(is_dna) != 1)) {
      stop("'is_dna' should be TRUE, FALSE or NULL")
    }
    if ((is_dna && has_U) ||
        (!is_dna && has_T)) {
      stop("if 'is_dna' is TRUE, sequences cannot contain 'U'; if is FALSE, sequences cannot contain 'T'")
    }
    if (!is_dna) dict["A"] <- "U"
  }
  
  inds_fun <- match(dict[alph], alph)
  names(inds_fun) <- as.character(1:length(alph))
  ret <- lapply(sq, function(s) pack_ints(inds_fun[s], alph_size))
  
  ret <- .set_alph(ret, alph)
  .set_class(ret, "nuc", TRUE)
}

