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

.bitify_sq <- function(sq, alph) {
  sq <- .char_to_int(sq, alph) 
  alph_size <- .get_alph_size(alph)
  lapply(sq, function(s) {
    .int_to_bit(s, alph_size)
  })
}

.bit_to_int <- function(s, alph_size) {
  if (length(s) == 1 && s == 0) {
    0L
  } else {
    s <- as.integer(unpack(s, alph_size))
    n <- length(s)
    tail_beg <- match(TRUE, (s == 0)) #this way of finding tail is sub-optimal, but works at least
    if (is.na(tail_beg)) {
      s
    } else {
      s[1:(tail_beg - 1)]
    }
  }
}

.debitify_sq <- function(sq, alph) {
  alph_size <- .get_alph_size(alph)
  na_val <- .get_na_val(alph)
  lapply(sq, function(s) {
    s <- .bit_to_int(s, alph_size)
    alph[s]
  })
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

.nc_bitify_sq_cnuc <- function(sq) {
  lapply(sq, function(s) {
    nc_pack_cnuc(charToRaw(s))
  })
}

.debitify_sq_cnuc <- function(sq) {
  na_char <- .get_na_char()
  sapply(sq, function(s) {
    rawToChar(unpack_nc_cnuc(s, na_char))
  })
}