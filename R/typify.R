#' Set type of untyped and atypical sequences
#' 
#' Function to set type of \strong{atp} or \strong{unt} sequences to \strong{ami} or 
#' \strong{nuc} provided that they do not have invalid letters (which can be checked
#' using \code{\link{get_sq_alphabet}} and \code{\link{get_invalid_letters}}).
#' 
#' @param sq an object of class \code{\link{sq}} with one of the types \strong{ami}, 
#' \strong{nuc}, \strong{unt} or \strong{atp}.
#' @param dest_type \code{\link{character}} string, destination type, either "ami" or "nuc".
#' 
#' @return An object of class \code{sq} that represents the same sequences as input \code{sq},
#' but with type as specified in \code{dest_type}.
#' 
#' @details 
#' Sometimes during reading \code{sq} from fasta file (using \code{\link{read_fasta}}) or 
#' constructing from character vector (using \code{\link{construct_sq}} or 
#' \code{\link{as.character}}) there are letters that are not in standard alphabets (see
#' \code{\link{aminoacids_df}} and \code{\link{nucleotides_df}}). In consequence, newly created
#' \code{sq} objects might have other type - \strong{atp} or \strong{unt}. After removal of those
#' non-standard letters (using \code{\link{substitute_letters}}), user might want to set type of 
#' \code{sq} object to one of standard types - \strong{ami} and \strong{nuc} - this is demanded
#' by some functions. This is what this function is designed for.
#' 
#' If \code{dest_type} is equal to type of \code{sq}, function does not do anything.
#' 
#' If \code{sq} object contain letters that are not in specified destination alphabet, an error 
#' will be thrown. However, \code{typify} converts automatically lowercase letters to uppercase
#' letters, so e.g. "a" is treated as element of amino acid alphabet.
#' 
#' If \code{sq} contains both lower and uppercase letters, they will be converted to uppercase, but
#' a message informing about it will be printed in the console. This action is default and can
#' be changed in package options (see \code{\link{tidysq-options}}).
#' 
#' Output \code{sq} object will not have \strong{cln} subtype, even if all letters of it fit in
#' clean alphabet of destination type (with exception of passing already clean object as input).
#' 
#' @examples 
#' # constructing artificial sq object with strange characters - type will be "unt"
#' sq_unt <- construct_sq(c("&VPLG&#", "##LCG"))
#' 
#' # substituting letters with "X", which stands for unknown amino acid
#' sq_sub <- substitute_letters(sq_unt, c(`&` = "X", `#` = "X"))
#' 
#' # setting type
#' typify(sq_sub, "ami")
#' 
#' @export 
typify <- function(sq, dest_type) {
  .validate_sq(sq)
  .check_type(dest_type, "'dest_type'")
  type <- .get_sq_type(sq)
  if (type == dest_type) {
    return(sq)
  }
  alph <- .get_alph(sq)
  up_alph <- unique(toupper(alph))
  dest_alph <- .get_standard_alph(dest_type, FALSE)
  
  .check_all_up_alph_proper(up_alph, dest_alph)
  if (!(length(alph) == length(up_alph))) {
    .handle_opt_txt("tidysq_a_typify_small_cap_let",
                    "in 'alphabet' attribute of 'sq' some letters appear as both lower and capital")
  }
  
  pack_fun <- switch(dest_type,
    "ami" = nc_pack_ami,
    "dna" = nc_pack_dna,
    "rna" = nc_pack_rna
  )
  ret <- .apply_sq(sq, "char", "char", function(s) {
    toupper(s)
  }, im_alph = dest_alph)
  
  ret <- .set_alph(ret, dest_alph)
  .set_class(ret, dest_type)
}
