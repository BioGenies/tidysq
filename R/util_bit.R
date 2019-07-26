.char_to_int <- function(sq, alph) {
  sq <- strsplit(sq, "")
  na_val <- .get_na_val(alph)
  sq <- lapply(sq, function(s) {
    s <- match(s, alph)
    if (length(s) == 0) s <- 0
    s[is.na(s)] <- na_val
    s
  })
  sq
}

.int_to_bit <- function(s, alph_size) {
  if (length(s) == 1 && s == 0) {
    as.raw(0)
  } else {
    ne <- length(s) / 8
    if (floor(ne) != ceiling(ne)) s[(length(s) + 1):(ceiling(ne) * 8)] <- 0
    pack(s, alph_size)
  }
}

.recode_sq <- function(sq, alph, new_alph, inds_func) {
  alph_size <- .get_alph_size(alph)
  new_alph_size <- .get_alph_size(new_alph)
  na_val <- .get_na_val(alph)
  new_na_val <- .get_na_val(new_alph)
  lapply(sq, function(s) {
    s <- .bit_to_int(s, alph_size)
    n <- length(s)
    s[s == na_val] <- NA
    tail_beg <- match(TRUE, (s == 0)) 
    if (!is.na(tail_beg)) {
      if (n == 1) {
        return(as.raw(0))
      } else {
        s <- s[1:(tail_beg-1)]
      }
    } 
    s <- inds_func[as.character(s)]  
    s[is.na(s)] <- new_na_val
    .int_to_bit(s, new_alph_size)
  })
}

.bitify_sq <- function(sq, type = NULL, is_clean = NULL, alph = NULL) {
  if (is.null(alph)) {
    pack_fun <- if (type == "ami") {
      if (is_clean) pack_cami else pack_ami
    } else if (type == "nuc") {
      if (is_clean) pack_cnuc else pack_nuc
    }
    lapply(sq, function(s) {
      pack_fun(charToRaw(s))
    })
  } else {
    alph <- charToRaw(paste(alph, collapse = ""))
    alph_size <- .get_alph_size(alph)
    lapply(sq, function(s) {
      pack_unt(charToRaw(s), alph, alph_size)
    })
  }
}

.debitify_sq <- function(sq) {
  na_char <- .get_na_char()
  type <- .get_sq_type(sq)
  if (type %in% c("ami", "nuc")) {
    is_clean <- .is_cleaned(sq)
    unpack_fun <- if (type == "ami") {
      if (is_clean) unpack_cami else unpack_ami
    } else if (type == "nuc") {
      if (is_clean) unpack_cnuc else unpack_nuc
    }
    sapply(sq, function(s) {
      rawToChar(unpack_fun(s, na_char))
    })
  } else {
    alph <- charToRaw(paste(.get_alph(sq), collapse = ""))
    alph_size <- .get_alph_size(alph)
    sapply(sq, function(s) {
      rawToChar(unpack_unt(s, na_char, alph, alph_size))
    })
  }
}