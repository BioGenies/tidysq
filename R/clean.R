#' Clean sequences containing ambiguous elements
#' 
#' Removes sequences containing ambiguous elements or removes ambiguous 
#' elements from sequences in a \code{sq} object.
#' 
#' @param sq a \code{\link{sq}} object of type \strong{ami}, \strong{dna} or \strong{rna}.
#' @param only_elements \code{\link{logical}} value indicating if only ambiguous elements
#' (i.e., matching more than one amino acid/nucleotide) of sequences should
#' be removed. If \code{FALSE} (default) whole sequences containing ambiguous 
#' elements are removed.
#'  
#' @return A \code{\link{sq}} object with the \strong{cln} subtype. 
#' 
#' @details This function allows cleaning of sequences containing ambiguous
#' elements. By default, sequences containing ambiguous elements are removed 
#' and \code{\link[=sq]{NULL}} (empty) sequences are introduced in their place. If 
#' \code{only_elements = TRUE} then only ambiguous elements are removed 
#' from sequences in \code{sq} object. Ambiguous letters of the amino acid alphabet 
#' are: B, J, O, U, X, Z, and these of the nucleotide alphabets are: W, S, M, K, R, 
#' Y, B, D, H, V, N. They are marked as 'amb' in \code{\link{aminoacids_df}} 
#' and \code{\link{nucleotides_df}} respectively. \code{\link[=sq]{NULL}} values (empty sequences) 
#' can be identified using \code{\link{is_null_sq}} function. 
#'
#' @examples 
#' # Creating objects to work on:
#' sq_ami <- construct_sq(c("MIAANYTWIL","TIAALGNIIYRAIE", 
#'                          "NYERTGHLI", "MAYXXXIALN"), type = "ami")
#' sq_dna <- construct_sq(c("ATGCAGGA", "GACCGAACGAN", 
#'                          "TGACGAGCTTA", "ACTNNAGCN"), type = "dna")
#' 
#' # Removing sequences with ambiguous elements:
#' clean(sq_ami)
#' clean(sq_dna)
#' 
#' # Removing ambiguous elements from sequences:
#' clean(sq_ami, only_elements = TRUE)
#' clean(sq_dna, only_elements = TRUE)
#' 
#' # Testing for presence of empty sequences after cleaning:
#' cln_sq <- clean(sq_ami)
#' is_null_sq(cln_sq)
#' 
#' 
#' @seealso \code{\link{sq}} \code{\link{aminoacids_df}} \code{\link{nucleotides_df}}
#' \code{\link{is_null_sq}}
#' @export
clean <- function(sq, only_elements = FALSE) {
  .validate_sq(sq)
  type <- .get_sq_type(sq)
  is_clean <- .is_cleaned(sq)
  .check_logical(only_elements, "'only_elements'", single_elem = TRUE)
  .check_type(type, "type of 'sq' object")
  
  if (is_clean) {
    return(sq)
  }
  alph_cln <- .get_standard_alph(type, TRUE)
  
  if (only_elements) {
    ret <- .apply_sq(sq, "char", "char", function(s) {
      s[(s %in% alph_cln)]
    }, alph_cln)
  } else {
    ret <- .apply_sq(sq, "char", "char", function(s) {
      if (!all(s %in% alph_cln)) raw(0) else s
    }, alph_cln) 
  }
  
  new_list_of(ret,
              ptype = raw(),
              alphabet = alph_cln,
              class = c(.get_sq_subclass(sq), "clnsq", "sq"))
}