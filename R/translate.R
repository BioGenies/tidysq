#' @export
translate <- function(sq, table = 1, ...) {
  UseMethod("translate")
}

#' @export
translate.default <- function(sq, table = 1, ...) {
  stop("cannot translate something that is neither DNA nor RNA sequence", call. = FALSE)
}

#' @export
translate.dnasq <- function(sq, table = 1, ...) {
  if (!.is_cleaned(sq)) {
    stop("sequence has to be cleaned to translate", call. = FALSE)
  }
  construct_sq_ami(Cpp_translate(as.character(sq), table),
                   is_clean = TRUE)
}

#' @export
translate.rnasq <- function(sq, table = 1, ...) {
  if (!.is_cleaned(sq)) {
    stop("sequence has to be cleaned to translate", call. = FALSE)
  }
  # a hack to avoid creating duplicate codon tables, at least for now
  # optimally should be deleted once the code operates without unpacking
  sq <- substitute_letters(sq, c(U = "T"))
  construct_sq_ami(Cpp_translate(as.character(sq), table),
                   is_clean = TRUE)
}
