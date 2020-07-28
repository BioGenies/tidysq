#' @export
vec_ptype_abbr.amisq <- function(x, ...) "ami"

#' @export
vec_ptype_abbr.dnasq <- function(x, ...) "dna"

#' @export
vec_ptype_abbr.rnasq <- function(x, ...) "rna"

#' @export
vec_ptype_abbr.untsq <- function(x, ...) "unt"

#' @export
vec_ptype_abbr.atpsq <- function(x, ...) "atp"

#' @export
vec_ptype_abbr.encsq <- function(x, ...) "enc"

#' @export
vec_ptype_abbr.clnsq <- function(x, ...) paste0("(c)", NextMethod())

#' @export
format.sq <- function(x, ...) {
  "hidden"
}
