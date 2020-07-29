#' @export
vec_cast.dnasq.dnasq <- function(x, to, ...) {
  if (.is_cleaned(x) == .is_cleaned(to)) x else construct_sq_dna(as.character(x), .is_cleaned(to))
}
#' @export
vec_cast.dnasq.character <- function(x, to, ...) construct_sq_dna(x, .is_cleaned(to))

#' @export
vec_cast.rnasq.rnasq <- function(x, to, ...) {
  if (.is_cleaned(x) == .is_cleaned(to)) x else construct_sq_rna(as.character(x), .is_cleaned(to))
}
#' @export
vec_cast.rnasq.character <- function(x, to, ...) construct_sq_rna(x, .is_cleaned(to))
