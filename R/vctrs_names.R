#' @export
vec_ptype_abbr.sq <- function(x, ...) ""

#' @export
vec_ptype_full.sq <- function(x, ...) ""

#' @export
vec_ptype_abbr.clnsq <- function(x, ...) "(c)"

#' @export
vec_ptype_full.clnsq <- function(x, ...) ", cln (cleaned)"

#' @export
vec_ptype_abbr.amisq <- function(x, ...) paste0(NextMethod(), "ami")

#' @export
vec_ptype_full.amisq <- function(x, ...) paste0("ami (amino acids)", NextMethod())

#' @export
vec_ptype_abbr.dnasq <- function(x, ...) paste0(NextMethod(), "dna")

#' @export
vec_ptype_full.dnasq <- function(x, ...) paste0("dna (DNA)", NextMethod())

#' @export
vec_ptype_abbr.rnasq <- function(x, ...) paste0(NextMethod(), "rna")

#' @export
vec_ptype_full.rnasq <- function(x, ...) paste0("rna (RNA)", NextMethod())

#' @export
vec_ptype_abbr.untsq <- function(x, ...) paste0(NextMethod(), "unt")

#' @export
vec_ptype_full.untsq <- function(x, ...) paste0("unt (unspecified type)", NextMethod())

#' @export
vec_ptype_abbr.atpsq <- function(x, ...) paste0(NextMethod(), "atp")

#' @export
vec_ptype_full.atpsq <- function(x, ...) paste0("atp (atypical alphabet)", NextMethod())

#' @export
vec_ptype_abbr.encsq <- function(x, ...) paste0(NextMethod(), "enc")

#' @export
vec_ptype_full.encsq <- function(x, ...) paste0("enc (encoded values)", NextMethod())

#' @export
vec_ptype_abbr.sq_alphabet <- function(x, ...) "sq_alph"

#' @export
vec_ptype_full.sq_alphabet <- function(x, ...) "tidysq alphabet"
