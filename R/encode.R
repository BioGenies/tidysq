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
    if (!("-" %in% names(encoding))) {
      encoding <- c(encoding, `-` = NA)
      .handle_opt_txt("tidysq_encode_nogap_action",
                      "'encoding' don't contain '-', assuming NA")
    }
  }
  
  alph <- .get_alph(sq)
  if (!all(alph %in% names(encoding))) {
    if (length(encoding) == length(alph)) {
      names(encoding) <- alph
      .handle_opt_txt("tidysq_encode_noname_action",
                      "'encoding' is unnamed, assuming values corresponds to letters in order")
    } else stop("'encoding' should be a named vector which names are superset of alphabet of 'sq'")
  }
  
  attr(sq, "alphabet") <- encoding[alph]
  class(sq) <- c("encsq", "sq")
  sq
}