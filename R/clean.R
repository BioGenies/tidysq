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
#' and \code{\link{nucleotides_df}} respectively. \code{NULL sq} values can
#' be identified using \code{\link{is_null_sq}} function. 
#'
#' @examples 
#' # Creating objects to work on:
#' sq_ami <- construct_sq(c("MIAANYTWIL","TIAALGNIIYRAIE", 
#'                          "NYERTGHLI", "MAYXXXIALN"), type = "ami")
#' sq_nuc <- construct_sq(c("ATGCAGGA", "GACCGAACGAN", 
#'                          "TGACGAGCTTA", "ACTNNAGCN"), type = "nuc")
#' 
#' # Removing sequences with ambiguous elements:
#' clean(sq_ami)
#' clean(sq_nuc)
#' 
#' # Removing ambiguous elements from sequences:
#' clean(sq_ami, only_elements = TRUE)
#' clean(sq_nuc, only_elements = TRUE)
#' 
#' # Testing for presence of empty sequences after cleaning:
#' cln_sq <- clean(sq_ami)
#' is_null_sq(cln_sq)
#' 
#' 
#' @seealso sq clnsq aminoacids_df nucleotides_df is_null_sq
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
  alph_cln <- .get_standard_alph(type, TRUE)
  
  
  if (only_elements) {
    ret <- .apply_sq(sq, "char", "none", function(s) {
      pack_chars(s[(s %in% alph_cln)], alph_cln)
    })
  } else {
    ret <- .apply_sq(sq, "char", "none", function(s) {
      if (!all(s %in% alph_cln)) raw(1) else
        pack_chars(s, alph_cln)
    }) 
  }

  class(ret) <- c("clnsq", class(sq))
  .set_alph(ret, alph_cln)
}