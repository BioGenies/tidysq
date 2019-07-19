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