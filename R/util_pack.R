pack <- function(x, alphabet,
                 NA_letter = getOption("tidysq_NA_letter"),
                 ignore_case = FALSE) {
  UseMethod("pack")
}

pack.character <- function(x, alphabet,
                           NA_letter = getOption("tidysq_NA_letter"),
                           ignore_case = FALSE) {
  CPP_pack_STRING(x, alphabet, NA_letter, ignore_case)
}

pack.list <- function(x, alphabet,
                      NA_letter = getOption("tidysq_NA_letter"),
                      ignore_case = FALSE) {
  pack_fun <- CPP_pack_RAWS
  if (length(x) > 0) {
    pack_fun <- switch (mode(x[[1]]),
                        character = CPP_pack_STRINGS,
                        numeric = CPP_pack_INTS,
                        raw = CPP_pack_RAWS
    )
  }
  pack_fun(x, alphabet, NA_letter, ignore_case)
}
