#' Create complement sequence from nucseq object 
#' 
#' @description Create complementary sequence from given nucleoide sequence 
#' according to complementarity.
#' 
#' @param nucsq \code{\link{nucsq}} object of type 'nuc'.
#' @param is_dna logical indicating if entered sequence is dna or rna. If 'is_dna' is \code{TRUE}, sequences cannot contain 'U' (uracil);
#' if is \code{FALSE}, sequences cannot contain 'T' (thymine). If \code{NULL} (default) the sequence type is not specified.
#'
#' @return \code{\link{nucsq}} object of the same type as input nucsq but built of complementary nucleotides to entered sequence.
#' 
#' @details This function allows obtaining  complement sequence which is created by 
#' matching elements (nucleotides) with complementarity to input nucsq object. If \code{is_dna = TRUE} entered sequence is dna. 
#' If \code{is_dna = FALSE} entered sequence is rna. By default the sequence type is not specified and function
#' needs to guess which type of sequence was entered. If sequence contain 'U' without 'T' the type will be set to rna.
#' If sequence contain 'T' without 'U' the type will be set to dna. If both 'U' and 'T' are present in the sequence
#' or sequence contains 'A' (adenine), but doesn't contain neither 'T' nor 'U'
#' or logical indicating is typed wrong (i.e., when sequence contain 'U' and logical indicating is set to dna) error will appear.
#' Both rna and dna sequences can be rewritten to complementary sequence. 
#' 
#' @examples 
#' # creating objects to work on:
#' nuc_dna_sequence <- construct_sq(c("ACTGCTG", "CTTAGA", "CCCT", "CTGAATGT"), type = "nuc")
#' nuc_rna_sequence <- construct_sq(c("ACUGCUG", "CUUAGA", "CCCU", "CUGAAUGU"), type = "nuc")
#' 
#' 
#' #creating complementary sequences with defined sequence type:
#' complement(nuc_dna_sequence, is_dna = TRUE)
#' complement(nuc_rna_sequence, is_dna = FALSE)
#' 
#' #creating complementary sequences without defined sequence type:
#' complement(nuc_dna_sequence)
#' complement(nuc_rna_sequence)
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
  sq <- lapply(nucsq, function(s) .bit_to_int(s, alph_size))
  
  has_U <- any(unlist(sq) == match("U", alph))
  has_T <- any(unlist(sq) == match("T", alph))
  has_A <- any(unlist(sq) == match("A", alph))
  
  if (has_U && has_T) {
    stop("'nucsq' sequences contains both 'U' and 'T' letters - should contain only one of them")
  }
  
  dict <- c(G = "C", C = "G", T = "A", U = "A", A = "T", `-` = "-")
  if (is.null(is_dna)) {
    if (!has_U && !has_T && has_A) {
      stop("'nucsq' sequences contains 'A' elements, but doesn't contain neither 'T' nor 'U' - unable to guess if it's dna or rna")
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
  
  inds_func <- match(dict[alph], alph)
  names(inds_func) <- as.character(1:length(alph))
  ret <- lapply(sq, function(s) .int_to_bit(inds_func[s], alph_size))
  
  class(ret) <- c("clnsq", "nucsq", "sq")
  attr(ret, "alphabet") <- alph
  ret
}
