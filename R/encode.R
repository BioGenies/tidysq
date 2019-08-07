#' @exportClass encsq
#' @export
encode <- function(sq, encoding) {
  validate_sq(sq)
  type <- .get_sq_type(sq)
  
  if (!is.numeric(encoding) ||
      is.null(names(encoding))) {
    stop("'encoding' should be a named numeric vector")
  }
  
  if (length(unique(names(encoding))) != length(encoding))
    stop("there are non-unique names in 'encoding' vector")
  
  if (type %in% c("ami", "nuc"))
    names(encoding) <- toupper(names(encoding))
  
  alph <- .get_alph(sq)
  alph_size <- .get_alph_size(alph)
  is_given <- alph %in% names(encoding)
  if (!all(is_given)) {
    ind <- (1:length(alph))[!is_given]
    for (s in sq) {
      if (any(unpack_ints(s, alph_size) %in% ind)) {
        .handle_opt_txt("tidysq_encode_no_given_action",
                        "there are letters in alphabet of 'sq' that appear in sequences, but weren't given in 'encoding' - assuming NA")
        break
      }
    }
    non_given <- rep(NA_real_, length(ind))
    names(non_given) <- alph[ind]
    encoding <- c(encoding, non_given)
  }
  
  sq <- .set_alph(sq, encoding[alph])
  .set_class(sq, "enc", FALSE)
}