# amisq ----
#' @export
vec_ptype2.amisq.amisq <- function(x, y, ...)
  .construct_sq_ptype("ami", is_clean = .is_cleaned(x) & .is_cleaned(y))
#' @export
vec_ptype2.amisq.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.amisq <- function(x, y, ...) y

# dnasq ----
#' @export
vec_ptype2.dnasq.dnasq <- function(x, y, ...)
  .construct_sq_ptype("dna", is_clean = .is_cleaned(x) & .is_cleaned(y))
#' @export
vec_ptype2.dnasq.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.dnasq <- function(x, y, ...) y

# rnasq ----
#' @export
vec_ptype2.rnasq.rnasq <- function(x, y, ...)
  .construct_sq_ptype("rna", is_clean = .is_cleaned(x) & .is_cleaned(y))
#' @export
vec_ptype2.rnasq.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.rnasq <- function(x, y, ...) y

# untsq ----
#' @export
vec_ptype2.untsq.untsq <- function(x, y, ...)
  .construct_sq_ptype("unt", alph = union(.get_alph(x), .get_alph(y)))

# atpsq ----
#' @export
vec_ptype2.atpsq.atpsq <- function(x, y, ...)
  .construct_sq_ptype("atp", alph = union(.get_alph(x), .get_alph(y)))
#' @export
vec_ptype2.atpsq.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.atpsq <- function(x, y, ...) y
