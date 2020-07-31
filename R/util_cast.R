# amisq ----
#' @export
vec_cast.amisq.amisq <- function(x, to, ...)
  if (.is_cleaned(x) == .is_cleaned(to)) x else construct_sq_ami(as.character(x), .is_cleaned(to))
#' @export
vec_cast.amisq.character <- function(x, to, ...) construct_sq_ami(x, .is_cleaned(to))
#' @export
vec_cast.character.amisq <- function(x, to, ...) unlist(.unpack_from_sq(x, "string"))

# dnasq ----
#' @export
vec_cast.dnasq.dnasq <- function(x, to, ...)
  if (.is_cleaned(x) == .is_cleaned(to)) x else construct_sq_dna(as.character(x), .is_cleaned(to))
#' @export
vec_cast.dnasq.character <- function(x, to, ...) construct_sq_dna(x, .is_cleaned(to))
#' @export
vec_cast.character.dnasq <- function(x, to, ...) unlist(.unpack_from_sq(x, "string"))

# rnasq ----
#' @export
vec_cast.rnasq.rnasq <- function(x, to, ...)
  if (.is_cleaned(x) == .is_cleaned(to)) x else construct_sq_rna(as.character(x), .is_cleaned(to))
#' @export
vec_cast.rnasq.character <- function(x, to, ...) construct_sq_rna(x, .is_cleaned(to))
#' @export
vec_cast.character.rnasq <- function(x, to, ...) unlist(.unpack_from_sq(x, "string"))

# untsq ----
#' @export
vec_cast.untsq.untsq <- function(x, to, ...)
#' @export
vec_cast.character.untsq <- function(x, to, ...) unlist(.unpack_from_sq(x, "string"))

# atpsq ----
#' @export
vec_cast.atpsq.atpsq <- function(x, to, ...)
  if (identical(.get_alph(x), .get_alph(to))) x else construct_sq(as.character(x), non_standard = .get_alph(to))
#' @export
vec_cast.atpsq.character <- function(x, to, ...) construct_sq(x, non_standard = .get_alph(to))
#' @export
vec_cast.character.atpsq <- function(x, to, ...) unlist(.unpack_from_sq(x, "string"))
