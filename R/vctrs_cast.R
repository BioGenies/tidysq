# amisq ----
#' @export
vec_cast.amisq.amisq <- function(x, to, ...)
  if (.is_cleaned(x) == .is_cleaned(to)) x else .construct_amisq(as.character(x), .is_cleaned(to))
#' @export
vec_cast.amisq.character <- function(x, to, ...) .construct_amisq(x, .is_cleaned(to))
#' @export
vec_cast.character.amisq <- function(x, to, ...) unlist(.unpack_from_sq(x, "string"))
#' @export
vec_cast.untsq.amisq <- function(x, to, ...) .construct_untsq(as.character(x), alph = alphabet(to))
#' @export
vec_cast.amisq.untsq <- function(x, to, ...) typify(x, "ami")

# dnasq ----
#' @export
vec_cast.dnasq.dnasq <- function(x, to, ...)
  if (.is_cleaned(x) == .is_cleaned(to)) x else .construct_dnasq(as.character(x), .is_cleaned(to))
#' @export
vec_cast.dnasq.character <- function(x, to, ...) .construct_dnasq(x, .is_cleaned(to))
#' @export
vec_cast.character.dnasq <- function(x, to, ...) unlist(.unpack_from_sq(x, "string"))
#' @export
vec_cast.untsq.dnasq <- function(x, to, ...) .construct_untsq(as.character(x), alph = alphabet(to))
#' @export
vec_cast.dnasq.untsq <- function(x, to, ...) typify(x, "dna")

# rnasq ----
#' @export
vec_cast.rnasq.rnasq <- function(x, to, ...)
  if (.is_cleaned(x) == .is_cleaned(to)) x else .construct_rnasq(as.character(x), .is_cleaned(to))
#' @export
vec_cast.rnasq.character <- function(x, to, ...) .construct_rnasq(x, .is_cleaned(to))
#' @export
vec_cast.character.rnasq <- function(x, to, ...) unlist(.unpack_from_sq(x, "string"))
#' @export
vec_cast.untsq.rnasq <- function(x, to, ...) .construct_untsq(as.character(x), alph = alphabet(to))
#' @export
vec_cast.rnasq.untsq <- function(x, to, ...) typify(x, "rna")

# untsq ----
#' @export
vec_cast.untsq.untsq <- function(x, to, ...)
  if (identical(alphabet(x), alphabet(to))) x else .construct_untsq(as.character(x), alph = alphabet(to))
#' @export
vec_cast.untsq.character <- function(x, to, ...) .construct_untsq(x, alph = alphabet(to))
#' @export
vec_cast.character.untsq <- function(x, to, ...) unlist(.unpack_from_sq(x, "string"))

# atpsq ----
#' @export
vec_cast.atpsq.atpsq <- function(x, to, ...)
  if (identical(alphabet(x), alphabet(to))) x else .nonst_construct_sq(as.character(x), alphabet(to))
#' @export
vec_cast.atpsq.character <- function(x, to, ...) .nonst_construct_sq(x, alphabet(to))
#' @export
vec_cast.character.atpsq <- function(x, to, ...) unlist(.unpack_from_sq(x, "string"))

# sq_alphabet ----
#' @export
vec_cast.sq_alphabet.character <- function(x, to, ...) sq_alphabet(x, na_character(to))
#' @export
vec_cast.character.sq_alphabet <- function(x, to, ...) vec_data(x)
