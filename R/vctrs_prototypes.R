# amisq ----
#' @export
vec_ptype2.amisq.amisq <- function(x, y, ...)
  .construct_sq_ptype("ami", is_clean = .is_cleaned(x) & .is_cleaned(y))
#' @export
vec_ptype2.amisq.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.amisq <- function(x, y, ...) y
#' @export
vec_ptype2.amisq.untsq <- function(x, y, ...)
  .construct_sq_ptype("unt", alph = union(alphabet(x), alphabet(y)))
#' @export
vec_ptype2.untsq.amisq <- function(x, y, ...)
  .construct_sq_ptype("unt", alph = union(alphabet(x), alphabet(y)))

# dnasq ----
#' @export
vec_ptype2.dnasq.dnasq <- function(x, y, ...)
  .construct_sq_ptype("dna", is_clean = .is_cleaned(x) & .is_cleaned(y))
#' @export
vec_ptype2.dnasq.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.dnasq <- function(x, y, ...) y
#' @export
vec_ptype2.dnasq.untsq <- function(x, y, ...)
  .construct_sq_ptype("unt", alph = union(alphabet(x), alphabet(y)))
#' @export
vec_ptype2.untsq.dnasq <- function(x, y, ...)
  .construct_sq_ptype("unt", alph = union(alphabet(x), alphabet(y)))

# rnasq ----
#' @export
vec_ptype2.rnasq.rnasq <- function(x, y, ...)
  .construct_sq_ptype("rna", is_clean = .is_cleaned(x) & .is_cleaned(y))
#' @export
vec_ptype2.rnasq.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.rnasq <- function(x, y, ...) y
#' @export
vec_ptype2.rnasq.untsq <- function(x, y, ...)
  .construct_sq_ptype("unt", alph = union(alphabet(x), alphabet(y)))
#' @export
vec_ptype2.untsq.rnasq <- function(x, y, ...)
  .construct_sq_ptype("unt", alph = union(alphabet(x), alphabet(y)))

# untsq ----
#' @export
vec_ptype2.untsq.untsq <- function(x, y, ...)
  .construct_sq_ptype("unt", alph = union(alphabet(x), alphabet(y)))
#' @export
vec_ptype2.untsq.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.untsq <- function(x, y, ...) y

# atpsq ----
#' @export
vec_ptype2.atpsq.atpsq <- function(x, y, ...)
  .construct_sq_ptype("atp", alph = union(alphabet(x), alphabet(y)))
#' @export
vec_ptype2.atpsq.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.atpsq <- function(x, y, ...) y

# sq_alphabet ----
#' @export
vec_ptype2.sq_alphabet.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.sq_alphabet <- function(x, y, ...) y
