.bitify_sq <- function(sq, alph) {
  if (is.numeric(sq[[1]])) 
    pack_fun <- pack_ints
  else if (any(lengths(sq) > 1)) 
    pack_fun <- pack_chars
  else pack_fun <- function(s, alph) pack_string(charToRaw(s), alph)
  
  alph_size <- .get_alph_size(alph)
  lapply(sq, function(s) {
    pack_fun(s, alph)
  })
}

.nc_bitify_sq <- function(sq, type, is_clean) {
  if      (type == "ami" &&  is_clean) pack_fun <- nc_pack_cami
  else if (type == "ami" && !is_clean) pack_fun <- nc_pack_ami
  else if (type == "nuc" &&  is_clean) pack_fun <- nc_pack_cnuc
  else if (type == "nuc" && !is_clean) pack_fun <- nc_pack_nuc
  
  lapply(sq, function(s) pack_fun(charToRaw(s)))
}

.debitify_sq <- function(sq, to) {
  alph <- .get_alph(sq)
  if (to == "char") 
    unpack_fun <- function(s) unpack_chars(s, alph, .get_na_char())
  else if (to == "int") 
    unpack_fun <- function(s) unpack_ints(s, .get_alph_size(alph))
  else if (to == "string") 
    unpack_fun <- function(s) unpack_string(s, alph, .get_na_char())
  
  lapply(sq, function(s) unpack_fun(s))
}

.apply_sq <- function(sq, ex_form, im_form, fun, im_alph = .get_alph(sq)) {
  ex_alph <- .get_alph(sq)
  if (ex_form == "char") 
    unpack_fun <- function(s) unpack_chars(s, ex_alph, .get_na_char())
  else if (ex_form == "int") 
    unpack_fun <- function(s) unpack_ints(s, .get_alph_size(ex_alph))
  else if (ex_form == "string") 
    unpack_fun <- function(s) unpack_string(s, ex_alph, .get_na_char())
 
  if (im_form == "char")
    pack_fun <- function(s) pack_chars(s, im_alph)
  else if (im_form == "int")
    pack_fun <- function(s) pack_ints(s, .get_alph_size(im_alph))
  else if (im_form == "string")
    pack_fun <- function(s) pack_string(charToRaw(s), im_alph)
  else if (im_form == "none")
    pack_fun <- identity

  lapply(sq, function(s) {
    s <- pack_fun(fun(unpack_fun(s)))
  })
}
