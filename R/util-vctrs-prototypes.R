# sq_dna_bsc ----
#' @export
vec_ptype2.sq_dna_bsc.sq_dna_bsc <- function(x, y, ...) x
#' @export
vec_ptype2.sq_dna_bsc.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.sq_dna_bsc <- function(x, y, ...) y

# sq_dna_ext ----
#' @export
vec_ptype2.sq_dna_ext.sq_dna_ext <- function(x, y, ...) x
#' @export
vec_ptype2.sq_dna_ext.sq_dna_bsc <- function(x, y, ...) x
#' @export
vec_ptype2.sq_dna_bsc.sq_dna_ext <- function(x, y, ...) y
#' @export
vec_ptype2.sq_dna_ext.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.sq_dna_ext <- function(x, y, ...) y

# sq_rna_bsc ----
#' @export
vec_ptype2.sq_rna_bsc.sq_rna_bsc <- function(x, y, ...) x
#' @export
vec_ptype2.sq_rna_bsc.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.sq_rna_bsc <- function(x, y, ...) y

# sq_rna_ext ----
#' @export
vec_ptype2.sq_rna_ext.sq_rna_ext <- function(x, y, ...) x
#' @export
vec_ptype2.sq_rna_ext.sq_rna_bsc <- function(x, y, ...) x
#' @export
vec_ptype2.sq_rna_bsc.sq_rna_ext <- function(x, y, ...) y
#' @export
vec_ptype2.sq_rna_ext.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.sq_rna_ext <- function(x, y, ...) y

# sq_ami_bsc ----
#' @export
vec_ptype2.sq_ami_bsc.sq_ami_bsc <- function(x, y, ...) x
#' @export
vec_ptype2.sq_ami_bsc.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.sq_ami_bsc <- function(x, y, ...) y

# sq_ami_ext ----
#' @export
vec_ptype2.sq_ami_ext.sq_ami_ext <- function(x, y, ...) x
#' @export
vec_ptype2.sq_ami_ext.sq_ami_bsc <- function(x, y, ...) x
#' @export
vec_ptype2.sq_ami_bsc.sq_ami_ext <- function(x, y, ...) y
#' @export
vec_ptype2.sq_ami_ext.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.sq_ami_ext <- function(x, y, ...) y

# sq_unt ----
#' @export
vec_ptype2.sq_unt.sq_unt <- function(x, y, ...)
  sq_ptype(union(as.character(alphabet(x)), as.character(alphabet(y))), "unt")
#' @export
vec_ptype2.sq_unt.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.sq_unt <- function(x, y, ...) y
#' @export
vec_ptype2.sq_unt.sq_dna_bsc <- function(x, y, ...)
  sq_ptype(union(as.character(alphabet(x)), as.character(alphabet(y))), "unt")
#' @export
vec_ptype2.sq_dna_bsc.sq_unt <- function(x, y, ...)
  sq_ptype(union(as.character(alphabet(x)), as.character(alphabet(y))), "unt")
#' @export
vec_ptype2.sq_unt.sq_dna_ext <- function(x, y, ...)
  sq_ptype(union(as.character(alphabet(x)), as.character(alphabet(y))), "unt")
#' @export
vec_ptype2.sq_dna_ext.sq_unt <- function(x, y, ...)
  sq_ptype(union(as.character(alphabet(x)), as.character(alphabet(y))), "unt")
#' @export
vec_ptype2.sq_unt.sq_rna_bsc <- function(x, y, ...)
  sq_ptype(union(as.character(alphabet(x)), as.character(alphabet(y))), "unt")
#' @export
vec_ptype2.sq_rna_bsc.sq_unt <- function(x, y, ...)
  sq_ptype(union(as.character(alphabet(x)), as.character(alphabet(y))), "unt")
#' @export
vec_ptype2.sq_unt.sq_rna_ext <- function(x, y, ...)
  sq_ptype(union(as.character(alphabet(x)), as.character(alphabet(y))), "unt")
#' @export
vec_ptype2.sq_rna_ext.sq_unt <- function(x, y, ...)
  sq_ptype(union(as.character(alphabet(x)), as.character(alphabet(y))), "unt")
#' @export
vec_ptype2.sq_unt.sq_ami_bsc <- function(x, y, ...)
  sq_ptype(union(as.character(alphabet(x)), as.character(alphabet(y))), "unt")
#' @export
vec_ptype2.sq_ami_bsc.sq_unt <- function(x, y, ...)
  sq_ptype(union(as.character(alphabet(x)), as.character(alphabet(y))), "unt")
#' @export
vec_ptype2.sq_unt.sq_ami_ext <- function(x, y, ...)
  sq_ptype(union(as.character(alphabet(x)), as.character(alphabet(y))), "unt")
#' @export
vec_ptype2.sq_ami_ext.sq_unt <- function(x, y, ...)
  sq_ptype(union(as.character(alphabet(x)), as.character(alphabet(y))), "unt")

# sq_atp ----
#' @export
vec_ptype2.sq_atp.sq_atp <- function(x, y, ...)
  sq_ptype(union(as.character(alphabet(x)), as.character(alphabet(y))), "atp")
#' @export
vec_ptype2.sq_atp.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.sq_atp <- function(x, y, ...) y

# sq_alphabet ----
#' @export
vec_ptype2.sq_alphabet.character <- function(x, y, ...) x
#' @export
vec_ptype2.character.sq_alphabet <- function(x, y, ...) y
