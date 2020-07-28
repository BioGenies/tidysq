#' @export
vec_ptype_abbr.amisq <- function(x, ...) "ami"

#' @export
vec_ptype_full.amisq <- function(x, ...) "ami (amino acids)"

#' @export
vec_ptype_abbr.dnasq <- function(x, ...) "dna"

#' @export
vec_ptype_full.dnasq <- function(x, ...) "dna (DNA)"

#' @export
vec_ptype_abbr.rnasq <- function(x, ...) "rna"

#' @export
vec_ptype_full.rnasq <- function(x, ...) "rna (RNA)"

#' @export
vec_ptype_abbr.untsq <- function(x, ...) "unt"

#' @export
vec_ptype_full.untsq <- function(x, ...) "unt (unspecified type)"

#' @export
vec_ptype_abbr.atpsq <- function(x, ...) "atp"

#' @export
vec_ptype_full.atpsq <- function(x, ...) "atp (atypical alphabet)"

#' @export
vec_ptype_abbr.encsq <- function(x, ...) "enc"

#' @export
vec_ptype_full.encsq <- function(x, ...) "enc (encoded values)"

#' @export
vec_ptype_abbr.clnsq <- function(x, ...) paste0("(c)", NextMethod())

#' @export
vec_ptype_full.clnsq <- function(x, ...) paste0("cln (cleaned), ", NextMethod())

#' @export
format.sq <- function(x, ...) {
  "hidden"
}
