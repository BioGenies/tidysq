#' Clean sequences containing ambiguous elements
#' 
#' Remove sequences containing ambiguous elements or remove ambiguous 
#' elements from sequences in a sq object.
#' 
#' @param sq \code{\link{sq}} object of type 'ami' or 'nuc'
#' @param only_elements logical indicating if only ambiguous elements
#' (i.e., matching more than one amino acid/nucleotide) of sequences should
#' be removed. If \code{FALSE} (default) whole sequences containing ambiguous 
#' elements are removed.
#'  
#' @return a \code{\link{sq}} object with the \code{\link{clnsq}} subtype. 
#' 
#' @details This function allows cleaning of sequences containing ambiguous
#' elements. By default, sequences containing ambiguous elements are removed 
#' and \code{NULL sq} values are introduced in their place. If 
#' \code{only_elements = TRUE} then only ambiguous elements are removed 
#' from sequences in sq object. Ambiguous letters of the amino acid alphabet 
#' are: B, J, O, U, X, Z, and of the nucleotide alphabet: W, S, M, K, R, 
#' Y, B, D, H, V, N. They are marked as 'amb' in \code{\link{aminoacids_df}} 
#' and \code{\link{nucleotides_df}} respectively.
#' 
#'
#' @examples 
#' # creating objects to work on:
#' sq_ami <- construct_sq(c("MIAANYTWIL","TIAALGNIIYRAIE", 
#'                          "NYERTGHLI", "MAYXXXIALN"), type = "ami")
#' sq_nuc <- construct_sq(c("ATGCAGGA", "GACCGAACGAN", 
#'                          "TGACGAGCTTA", "ACTNNAGCN"), type = "nuc")
#' 
#' # removing sequences with ambiguous elements:
#' clean(sq_ami)
#' clean(sq_nuc)
#' 
#' # removing ambiguous elements from sequences:
#' clean(sq_ami, only_elements = TRUE)
#' clean(sq_nuc, only_elements = TRUE)
#' 
#' @seealso sq clnsq aminoacids_df nucleotides_df
#' @exportClass clnsq
#' @export
clean <- function(sq, only_elements = FALSE) {
  validate_sq(sq)
  type <- .get_sq_type(sq)
  is_clean <- .is_cleaned(sq)

  if (!is.logical(only_elements) ||
      is.na(only_elements) ||
      length(only_elements) != 1) {
    stop("'only_elements' has to be either TRUE or FALSE")
  }
  if (!(type %in% c("ami", "nuc"))) {
    stop("function 'clean' is meant to be used only with 'ami' or 'nuc' sequences")
  }
  if (is_clean) {
    return(sq)
  }
  
  alph <- .get_alph(sq)
  alph_size <- .get_alph_size(alph)
  alph_cln <- if (type == "ami") 
    aminoacids_df[!aminoacids_df[["amb"]], "one"] 
  else
    nucleotides_df[!nucleotides_df[["amb"]], "one"]
  
  inds_uncln <- (1:length(alph))[!(alph %in% alph_cln)]
  
  if (only_elements) {
    ret <- lapply(sq, function(s) {
      s <- .bit_to_int(s, alph_size)
      .int_to_bit(s[!(s %in% inds_uncln)], alph_size)
    })
  } else {
    ret <- lapply(sq, function(s) {
      st <- .bit_to_int(s, alph_size)
      if (any(st %in% inds_uncln)) raw(1) else s
    })
  }
  
  inds_func <- 1:length(alph_cln)
  names(inds_func) <- as.character((1:length(alph))[-inds_uncln])
  
  ret <- .recode_sq(ret, alph, alph_cln, inds_func)
  
  class(ret) <- c("clnsq", class(sq))
  attr(ret, "alphabet") <- alph_cln
  ret
}