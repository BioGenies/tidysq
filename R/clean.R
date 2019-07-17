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