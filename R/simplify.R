#' @exportClass simsq
#' @export
simplify <- function(sq, encoding) {
  validate_sq(sq)
  type <- .get_sq_type(sq)
  
  if (!is.character(encoding) ||
      any(!(encoding %in% c(letters, "-")))) {
    stop("'encoding' should be a named vector, where names are upper latin letters or '-' and elements are lower latin letters and '-'")
  }
  
  alph <- .get_alph(sq)
  if (!all(alph %in% names(encoding))) {
    stop("each of letters (elements of 'alphabet' attribute of 'sq') should appear in names of encoding")
  }
  
  new_alph <- sort(unique(encoding))
  inds_func <- match(encoding, new_alph)
  names(inds_func) <- as.character(match(names(encoding), alph))
  
  ret <- lapply(sq, function(s) inds_func[s])
  class(ret) <- c("simsq", "sq")
  attr(ret, "alphabet") <- new_alph
  ret
}

