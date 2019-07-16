.get_alph_size <- function(alph) {
  2 ^ ceiling(log(length(alph) + 2, 4))
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
  n <- length(s)
  if (alph_size == 8 || 
      n < 2) {
    as.raw(s)
  } else {
    #s <- s * rep(2 ^ (0:(8 / alph_size - 1) * alph_size), length.out = n)
    
    if (alph_size == 4) {
      if (n %% 2 == 1) {
        s <- c(s, 0)
        n <- n + 1
      }
      as.raw(s[seq(1, n, by = 2)] + 16 * s[seq(2, n, by = 2)])
    } else if (alph_size == 8 &&
               n < 4) {
      as.raw(sum(s))
    } else {
      # WIP
      # as.raw(s[seq(1, n, by = 2)] + 
      #          s[seq(2, n, by = 2)] +
      #          s[seq(3, n, by = 2)] +
      #          s[seq(4, n, by = 2)])
    }
    
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
  if (alph_size == 8) {
    as.integer(s)
  } else if (alph_size == 4) {
    n <- 2 * length(s)
    ret <- integer(n)
    s <- as.integer(s)
    ret[seq(1, n, by = 2)] <- s %% 16
    ret[seq(2, n, by = 2)] <- s %/% 16
    ret
  } else if (alph_size == 2) {
    n <- 4 * length(s)
    ret <- integer(n)
    s <- as.integer(s)
    ret[seq(1, n, by = 4)] <- s %% 4
    ret[seq(2, n, by = 4)] <- (s %/% 4) %% 4
    ret[seq(3, n, by = 4)] <- (s %/% 16) %% 4
    ret[seq(4, n, by = 4)] <- s %/% 64
    ret
  }
}

.debitify_sq <- function(sq, alph) {
  alph_size <- .get_alph_size(alph)
  na_val <- .get_na_val(alph)
  lapply(sq, function(s) {
    s <- .bit_to_int(s, alph_size)
    n <- length(s)
    s[s == na_val] <- NA
    tail_beg <- match(TRUE, (s[(n - 8 / alph_size + 1):n] == 0)) 
    if (is.na(tail_beg)) {
      alph[s]
    } else {
      alph[s[1:(n - 8 / alph_size + tail_beg - 1)]]
    }
  })
}
