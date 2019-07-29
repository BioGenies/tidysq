.bitify_sq <- function(sq, type = NULL, is_clean = NULL, alph = NULL) {
  if (is.numeric(sq[[1]])) 
    t_fun <- rawize
  else if (any(lengths(sq) > 1)) 
    t_fun <- function(x) charToRaw(paste(x, collapse = ""))
  else t_fun <- charToRaw
  if (is.null(alph)) {
    pack_fun <- if (type == "ami") {
      if (is_clean) pack_cami else pack_ami
    } else if (type == "nuc") {
      if (is_clean) pack_cnuc else pack_nuc
    }
    lapply(sq, function(s) {
      pack_fun(t_fun(s))
    })
  } else {
    alph_size <- .get_alph_size(alph)
    alph <- .rawize_alph(alph)
    lapply(sq, function(s) {
      pack_unt(t_fun(s), alph, alph_size)
    })
  }
}

.debitify_sq <- function(sq, to = "char") {
  if (to == "char") 
    t_fun <- rawToChar
  else if (to == "int") 
    t_fun <- integerize
  else if (to == "chars") 
    t_fun <- function(x) strsplit(rawToChar(x), "")[[1]]
  
  na_char <- .get_na_char()
  type <- .get_sq_type(sq)
  if (type %in% c("ami", "nuc")) {
    is_clean <- .is_cleaned(sq)
    unpack_fun <- if (type == "ami") {
      if (is_clean) unpack_cami else unpack_ami
    } else if (type == "nuc") {
      if (is_clean) unpack_cnuc else unpack_nuc
    }
    sapply(sq, function(s) t_fun(unpack_fun(s, na_char)))
  } else {
    alph <- charToRaw(paste(.get_alph(sq), collapse = ""))
    alph_size <- .get_alph_size(alph)
    sapply(sq, function(s) t_fun(unpack_unt(s, na_char, alph, alph_size)))
  }
}


.recode_sq <- function(sq, alph, new_alph, letters_fun) {
  na_val <- .get_na_val(alph)
  new_na_val <- .get_na_val(new_alph)
  alph_size <- .get_alph_size(alph)
  new_alph_size <- .get_alph_size(new_alph)
  alph <- .rawize_alph(alph)
  new_alph <- .rawize_alph(new_alph)
  lapply(sq, function(s) {
    s <- repack(s, na_char, alph, 
                new_alph, 
                alph_size, new_alph_size, 
                letters_fun)
  })
}

.apply_sq <- function(sq, intern_form, fun) {
  if (intern_form == "char") {
    t_fun <- rawToChar
  } else if (intern_form == "int") {
    t_fun <- integerize
  } else if (intern_form == "chars") {
    t_fun <- function(x) strsplit(rawToChar(x), "")[[1]]
  } 
  alph <- .get_alph(sq)
  alph_size <- .get_alph_size(alph)
  na_char <- .get_na_char()
  alph <- .rawize_alph(alph)
  lapply(sq, function(s) {
    s <- t_fun(unpack_unt(s, na_char, alph, alph_size))
    fun(s)
  })
}
