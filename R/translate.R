#' @export
translate <- function(sq, ...) {
  UseMethod("translate")
}

#' @export
translate.default <- function(sq, ...) {
  stop("cannot translate something that is neither DNA nor RNA sequence", call. = FALSE)
}

#' @export
translate.dnasq <- function(sq, ...) {
  if (!.is_cleaned(sq)) {
    stop("sequence has to be cleaned to translate", call. = FALSE)
  }
  sq <- as.character(sq)
  codons <- lapply(sq, extractCodons)
  ret <- vapply(codons, codonsToAminoAcids, character(1))
  construct_sq_ami(ret, is_clean = TRUE)
}
