.get_alph_size <- function(alph) {
  ceiling(log2(length(alph) + 2))
}

.get_na_val <- function(alph) {
  2 ^ .get_alph_size(alph) - 1
}

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
    pack(s, alph_size)
  }
}

.int_to_bit2 <- function(s, alph_size) {
  if (length(s) == 1 && s == 0) {
    as.raw(0)
  } else {
    ne <- length(s) / 8
    eights <- lapply(1:floor(ne), function(ind) (8 * (ind - 1)):(8 * ind - 1))
    if (floor(ne) != ceiling(ne)) eights[[ceiling(ne)]] <- (floor(ne) * 8):length(s) 
    do.call(c, lapply(0:(eights - 1), function(i) pack(s[(8 * i):(8 * i + 7)], alph_size)))
  }
}

.bitify_sq <- function(sq, alph) {
  sq <- .char_to_int(sq, alph) 
  alph_size <- .get_alph_size(alph)
  lapply(sq, function(s) {
    .int_to_bit(s, alph_size)
  })
}

.bitify_sq <- function(sq, alph) {
  sq <- .char_to_int(sq, alph) 
  alph_size <- .get_alph_size(alph)
  lapply(sq, function(s) {
    .int_to_bit2(s, alph_size)
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
