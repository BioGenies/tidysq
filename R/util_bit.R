.pack_to_sq <- function(proto, alph) {
  if (length(proto) == 0) 
    list()
  else if (is.numeric(proto[[1]]))
    CPP_pack_INTS(proto, alph)
  else if (any(lengths(proto) > 1)) 
    CPP_pack_STRINGS(proto, alph)
  else 
    CPP_pack_STRING(proto, alph)
}

.nc_pack_to_sq <- function(proto, type, is_clean) {
  .pack_to_sq(proto, .get_standard_alph(type, is_clean))
}

.unpack_from_sq <- function(sq, to) {
  if (length(sq) == 0) 
    sq
  else if (to == "char") 
    CPP_unpack_STRINGS(sq)
  else if (to == "int") 
    CPP_unpack_INTS(sq)
  else if (to == "string") 
    CPP_unpack_STRING(sq)
}

.apply_sq <- function(sq, ex_form, im_form, fun, im_alph = alphabet(sq)) {
  if (length(sq) == 0) 
    return(sq)
  
  unpacking_fun <- switch (ex_form,
                           char = CPP_unpack_STRINGS,
                           int = CPP_unpack_INTS,
                           string = CPP_unpack_STRING
  )
  
  if (im_form == "none") 
    return(lapply(sq, function(s) fun(unpacking_fun(s))))
  
  packing_fun <- switch (im_form,
                         char = CPP_pack_STRINGS,
                         int = CPP_pack_INTS,
                         string = CPP_pack_STRING
  )

  packing_fun(lapply(1:length(sq), function(i) fun(unpacking_fun(sq[i])[[1]])), im_alph)
}
