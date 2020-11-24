# sq_dna_bsc ----
#' @export
vec_cast.sq_dna_bsc.sq_dna_bsc <- function(x, to, ...) x
#' @export
vec_cast.sq_dna_bsc.character <- function(x, to, ...)
  sq(x, alphabet = "dna_bsc")
#' @export
vec_cast.character.sq_dna_bsc <- function(x, to, ...)
  unpack(x, "STRING")
#' @export
vec_cast.sq_dna_bsc.sq_unt <- function(x, to, ...)
  typify(x, "dna_bsc")
#' @export
vec_cast.sq_unt.sq_dna_bsc <- function(x, to, ...)
  pack(as.character(x), alphabet(to))

# sq_dna_ext ----
#' @export
vec_cast.sq_dna_ext.sq_dna_ext <- function(x, to, ...) x
#' @export
vec_cast.sq_dna_ext.sq_dna_bsc <- function(x, to, ...)
  typify(x, "dna_ext")
#' @export
vec_cast.sq_dna_ext.character <- function(x, to, ...)
  sq(x, alphabet = "dna_ext")
#' @export
vec_cast.character.sq_dna_ext <- function(x, to, ...)
  unpack(x, "STRING")
#' @export
vec_cast.sq_dna_ext.sq_unt <- function(x, to, ...)
  typify(x, "dna_ext")
#' @export
vec_cast.sq_unt.sq_dna_ext <- function(x, to, ...)
  pack(as.character(x), alphabet(to))

# sq_rna_bsc ----
#' @export
vec_cast.sq_rna_bsc.sq_rna_bsc <- function(x, to, ...) x
#' @export
vec_cast.sq_rna_bsc.character <- function(x, to, ...)
  sq(x, alphabet = "rna_bsc")
#' @export
vec_cast.character.sq_rna_bsc <- function(x, to, ...)
  unpack(x, "STRING")
#' @export
vec_cast.sq_rna_bsc.sq_unt <- function(x, to, ...)
  typify(x, "rna_bsc")
#' @export
vec_cast.sq_unt.sq_rna_bsc <- function(x, to, ...)
  pack(as.character(x), alphabet(to))

# sq_rna_ext ----
#' @export
vec_cast.sq_rna_ext.sq_rna_ext <- function(x, to, ...) x
#' @export
vec_cast.sq_rna_ext.sq_rna_bsc <- function(x, to, ...)
  typify(x, "rna_ext")
#' @export
vec_cast.sq_rna_ext.character <- function(x, to, ...)
  sq(x, alphabet = "rna_ext")
#' @export
vec_cast.character.sq_rna_ext <- function(x, to, ...)
  unpack(x, "STRING")
#' @export
vec_cast.sq_rna_ext.sq_unt <- function(x, to, ...)
  typify(x, "rna_ext")
#' @export
vec_cast.sq_unt.sq_rna_ext <- function(x, to, ...)
  pack(as.character(x), alphabet(to))

# sq_ami_bsc ----
#' @export
vec_cast.sq_ami_bsc.sq_ami_bsc <- function(x, to, ...) x
#' @export
vec_cast.sq_ami_bsc.character <- function(x, to, ...)
  sq(x, alphabet = "ami_bsc")
#' @export
vec_cast.character.sq_ami_bsc <- function(x, to, ...)
  unpack(x, "STRING")
#' @export
vec_cast.sq_ami_bsc.sq_unt <- function(x, to, ...)
  typify(x, "ami_bsc")
#' @export
vec_cast.sq_unt.sq_ami_bsc <- function(x, to, ...)
  pack(as.character(x), alphabet(to))

# sq_ami_ext ----
#' @export
vec_cast.sq_ami_ext.sq_ami_ext <- function(x, to, ...) x
#' @export
vec_cast.sq_ami_ext.sq_ami_bsc <- function(x, to, ...)
  typify(x, "ami_ext")
#' @export
vec_cast.sq_ami_ext.character <- function(x, to, ...)
  sq(x, alphabet = "ami_ext")
#' @export
vec_cast.character.sq_ami_ext <- function(x, to, ...)
  unpack(x, "STRING")
#' @export
vec_cast.sq_ami_ext.sq_unt <- function(x, to, ...)
  typify(x, "ami_ext")
#' @export
vec_cast.sq_unt.sq_ami_ext <- function(x, to, ...)
  pack(as.character(x), alphabet(to))

# sq_unt ----
#' @export
vec_cast.sq_unt.sq_unt <- function(x, to, ...)
  if (identical(alphabet(x), alphabet(to))) x else pack(as.character(x), alphabet(to))
#' @export
vec_cast.sq_unt.character <- function(x, to, ...)
  pack(as.character(x), alphabet(to))
#' @export
vec_cast.character.sq_unt <- function(x, to, ...)
  unpack(x, "STRING")

# sq_atp ----
#' @export
vec_cast.sq_atp.sq_atp <- function(x, to, ...)
  if (identical(alphabet(x), alphabet(to))) x else pack(as.character(x), alphabet(to))
#' @export
vec_cast.sq_atp.character <- function(x, to, ...)
  pack(as.character(x), alphabet(to))
#' @export
vec_cast.character.sq_atp <- function(x, to, ...)
  unpack(x, "STRING")

# sq_alphabet ----
#' @export
vec_cast.sq_alphabet.sq_alphabet <- function(x, to, ...) {
  attr(x, "type") <- attr(to, "type")
  x
}
#' @export
vec_cast.sq_alphabet.character <- function(x, to, ...)
  sq_alphabet(x, attr(to, "type"))
#' @export
vec_cast.character.sq_alphabet <- function(x, to, ...)
  vec_data(x)
