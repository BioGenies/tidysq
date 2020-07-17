.pack_to_sq <- function(sq, alph) {
  if (length(sq) == 0) return(sq)
  if (is.numeric(sq[[1]])) 
    packing_fun <-  C_pack_ints
  else if (any(lengths(sq) > 1)) 
    packing_fun <- C_pack_chars
  else
    packing_fun <- function(s, alph) C_pack_string(charToRaw(s), alph)
  
  alph_size <- .get_alph_size(alph)
  ret <- lapply(sq, function(s) packing_fun(s, alph))
  .set_original_length(ret, nchar(sq))
}

.nc_pack_to_sq <- function(sq, type, is_clean) {
  if      (type == "ami" &&  is_clean) packing_fun <- nc_pack_cami
  else if (type == "ami" && !is_clean) packing_fun <- nc_pack_ami
  else if (type == "dna" &&  is_clean) packing_fun <- nc_pack_cdna
  else if (type == "dna" && !is_clean) packing_fun <- nc_pack_dna
  else if (type == "rna" &&  is_clean) packing_fun <- nc_pack_crna
  else if (type == "rna" && !is_clean) packing_fun <- nc_pack_rna
  
  ret <- lapply(sq, function(s) packing_fun(charToRaw(s)))
  .set_original_length(ret, nchar(sq))
}

.unpack_from_sq <- function(sq, to) {
  if (length(sq) == 0) return(sq)
  alph <- .get_alph(sq)
  if (to == "char") 
    unpacking_fun <- function(s) C_unpack_chars(s, alph, .get_na_char())
  else if (to == "int") 
    unpacking_fun <- function(s) C_unpack_ints(s, .get_alph_size(alph))
  else if (to == "string") 
    unpacking_fun <- function(s) C_unpack_string(s, alph, .get_na_char())
  
  ret <- lapply(sq, function(s) unpacking_fun(s))
  .set_original_length(ret, sapply(sq, attr, "original_length"))
}

.apply_sq <- function(sq, ex_form, im_form, fun, im_alph = .get_alph(sq)) {
  if (length(sq) == 0) return(sq)
  ex_alph <- .get_alph(sq)
  if (ex_form == "char") 
    unpacking_fun <- function(s) C_unpack_chars(s, ex_alph, .get_na_char())
  else if (ex_form == "int") 
    unpacking_fun <- function(s) C_unpack_ints(s, .get_alph_size(ex_alph))
  else if (ex_form == "string") 
    unpacking_fun <- function(s) C_unpack_string(s, ex_alph, .get_na_char())
  
  if (im_form == "none") return(lapply(sq, function(s) fun(unpacking_fun(s))))
  
  if (im_form == "char")
    packing_fun <- function(s) C_pack_chars(s, im_alph)
  else if (im_form == "int")
    packing_fun <- function(s) C_pack_ints(s, .get_alph_size(im_alph))
  else if (im_form == "string")
    packing_fun <- function(s) C_pack_string(charToRaw(s), im_alph)

  ret <- lapply(sq, function(s) fun(unpacking_fun(s)))
  orig_length <- lengths(ret)
  ret <- lapply(ret, function(s) packing_fun(s))
  .set_original_length(ret, orig_length)
}
