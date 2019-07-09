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
    stop("function 'clean' is meant to be used only with 'ami' or 'unt' sequences")
  }
  if (is_clean) {
    return(sq)
  }
  
  alph <- .get_alph(sq)
  alph_cln <- if (type == "ami") 
    aminoacids_df[!aminoacids_df[["amb"]], "one"] 
  else
    nucleotides_df[!nucleotides_df[["amb"]], "one"]
  
  inds_uncln <- (1:length(alph))[!(alph %in% alph_cln)]
  
  if (only_elements) {
    ret <- lapply(sq, function(s) s[!(s %in% inds_uncln)])
  } else {
    ret <- lapply(sq, function(s) if (any(s %in% inds_uncln)) integer(0) else s)
  }
  
  ret <- lapply(ret, function(s) match(alph[s], alph_cln))
  
  class(ret) <- c("clnsq", class(sq))
  attr(ret, "alphabet") <- alph_cln
  ret
}