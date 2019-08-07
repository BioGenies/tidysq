#' @exportClass encsq
#' @export
encode <- function(sq, encoding) {
  validate_sq(sq)
  type <- .get_sq_type(sq)
  
  if (!is.numeric(encoding)) {
    stop("'encoding' should be a numeric vector")
  }
  
  if (!is.null(names(encoding)) &&
      length(unique(names(encoding))) != length(encoding))
    stop("there are non-unique names in 'encoding' vector")
  
  if (type == "ami") {
    names(encoding) <- toupper(names(encoding))
    has_gap <- "-" %in% names(encoding)
    has_star <- "*" %in% names(encoding)
    if (!has_gap || !has_star) 
      .handle_opt_txt("tidysq_encode_nogap_action",
                      "'encoding' doesn't contain values for '-' or '*' or both, assuming NA")
    if (!has_gap) encoding <- c(encoding, `-` = NA)
    if (!has_star) encoding <- c(encoding, `*` = NA)
  } 
  if (type == "nuc") {
    if (!is.null(names(encoding))) names(encoding) <- toupper(names(encoding))
    if (!is.null(names(encoding))) {
      has_U <- "U" %in% names(encoding)
      has_T <- "T" %in% names(encoding)
      if (has_U && !has_T) encoding <- c(encoding, T = unname(encoding["U"]))
      if (has_T && !has_U) encoding <- c(encoding, U = unname(encoding["T"]))
    }
    if (.is_cleaned(sq) && is.null(names(encoding))) {
      if (length(encoding) == 4) encoding <- c(encoding, encoding[4])
      if (length(encoding) == 5) encoding <- c(encoding[c(1:4, 4)], `-` = encoding[5])
    } 
    has_gap <- "-" %in% names(encoding)
    if (!has_gap) {
      .handle_opt_txt("tidysq_encode_nogap_action",
                      "'encoding' doesn't contain values for '-', assuming NA")
      encoding <- c(encoding, `-` = NA)
    }
  }
  
  
  
  alph <- .get_alph(sq)
  if (!all(alph %in% names(encoding))) {
    if (length(encoding) == length(alph)) {
      names(encoding) <- alph
      .handle_opt_txt("tidysq_encode_noname_action",
                      "'encoding' is unnamed, assuming values corresponds to letters in order")
    } else stop("'encoding' should be a named vector which names are superset of alphabet of 'sq' or vector of appropriate length (see documentation for details)")
  }
  
  sq <- .set_alph(sq, encoding[alph])
  .set_class(sq, "enc", FALSE)
}