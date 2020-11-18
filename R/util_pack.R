pack <- function(x, alphabet,
                 NA_letter = getOption("tidysq_NA_letter"),
                 safe_mode = getOption("tidysq_safe_mode"),
                 ignore_case = FALSE) {
  UseMethod("pack")
}

pack.character <- function(x, alphabet,
                           NA_letter = getOption("tidysq_NA_letter"),
                           safe_mode = getOption("tidysq_safe_mode"),
                           ignore_case = FALSE) {
  # TODO: implement safe mode (exactly here)
  # TODO: passing x with letters that are not in alphabet throw an error
  #  e.g. passing "O" with "ami_bsc" type
  CPP_pack_STRING(x, alphabet, NA_letter, ignore_case)
}

pack.list <- function(x, alphabet,
                      NA_letter = getOption("tidysq_NA_letter"),
                      safe_mode = getOption("tidysq_safe_mode"),
                      ignore_case = FALSE) {
  # TODO: anyNA won't work with vctrs_list_of
  if (safe_mode && anyNA(x, recursive = TRUE))
    stop("NA value encountered during packing", call. = FALSE)
  # To avoid repeating the same args many times
  op <- CPP_pack_RAWS
  if (length(x) > 0) {
    op <- switch (mode(x[[1]]),
                 character = CPP_pack_STRINGS,
                 numeric = CPP_pack_INTS,
                 raw = CPP_pack_RAWS
    )
  }
  op(x, alphabet, NA_letter, ignore_case)
}
