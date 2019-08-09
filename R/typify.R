#' @export
typify <- function(sq, dest_type) {
  validate_sq(sq)
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
    .handle_opt_txt("tidysq_typify_small_cap_let",
                    "in 'alphabet' attribute of 'sq' some letters appear as both lower and capital")
  }
  
  pack_fun <- if (dest_alph == "ami") nc_pack_ami else nc_pack_nuc
  ret <- .apply_sq(sq, "char", "none", function(s) {
    pack_fun(s)
  })
  
  ret <- .set_alph(ret, dest_alph)
  .set_class(ret, dest_type)
}