#' @export
vec_ptype2.dnasq.dnasq <- function(x, y, ...) construct_sq_dna(character(), .is_cleaned(x) & .is_cleaned(y))
#' @export
vec_ptype2.dnasq.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.dnasq <- function(x, y, ...) y

#' @export
vec_ptype2.rnasq.rnasq <- function(x, y, ...) construct_sq_rna(character(), .is_cleaned(x) & .is_cleaned(y))
#' @export
vec_ptype2.rnasq.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.rnasq <- function(x, y, ...) y
