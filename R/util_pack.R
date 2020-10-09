pack <- function(x, alphabet, NA_letter, safe_mode) {
  UseMethod("pack")
}

pack.character <- function(x, alphabet, NA_letter, safe_mode) {
  # TODO: implement safe mode (exactly here)
  CPP_pack_STRING(x, alphabet, NA_letter)
}

pack.list <- function(x, alphabet, NA_letter, safe_mode) {
  # TODO: anyNA won't work with vctrs_list_of
  if (safe_mode && anyNA(x, recursive = TRUE))
    stop("NA value encountered during packing", call. = FALSE)
  # To avoid repeating the same args many times
  f <- CPP_pack_RAWS
  if (length(x) > 0) {
    f <- switch (mode(x[[1]]),
                 character = CPP_pack_STRINGS,
                 numeric = CPP_pack_INTS,
                 raw = CPP_pack_RAWS
    )
  }
  f(x, alphabet, NA_letter)
}
