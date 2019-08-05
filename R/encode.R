#' @exportClass encsq
#' @export
encode <- function(sq, encoding) {
  validate_sq(sq)
  type <- .get_sq_type(sq)
  
  if (!is.numeric(encoding)) {
    stop("'encoding' should be a numeric vector")
  }
  
  if (type %in% c("ami", "nuc")) {
    names(encoding) <- toupper(names(encoding))
    has_gap <- "-" %in% names(encoding)
    has_end <- "*" %in% names(encoding)
    
    if (!(has_gap && has_end)) {
      .handle_opt_txt("tidysq_encode_nogap_action",
                      "'encoding' don't contain values for '-' or '*' or both, assuming NA")
    }
    if (!has_gap) encoding <- c(encoding, `-` = NA)
    if (!has_end) encoding <- c(encoding, `*` = NA)
  }
  
  alph <- .get_alph(sq)
  if (!all(alph %in% names(encoding))) {
    if (length(encoding) == length(alph)) {
      names(encoding) <- alph
      .handle_opt_txt("tidysq_encode_noname_action",
                      "'encoding' is unnamed, assuming values corresponds to letters in order")
    } else stop("'encoding' should be a named vector which names are superset of alphabet of 'sq'")
  }
  
  sq <- .set_alph(sq, encoding[alph])
  .set_class(sq, "enc", FALSE)
}