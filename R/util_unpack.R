unpack <- function(x, format,
                   NA_letter = getOption("tidysq_NA_letter")) {
  assert_choice(format, c("RAWS", "INTS", "STRINGS", "STRING"))
  assert_string(NA_letter)
  
  op <- switch (
    format,
    RAWS = CPP_unpack_RAWS,
    INTS = CPP_unpack_INTS,
    STRINGS = CPP_unpack_STRINGS,
    STRING = CPP_unpack_STRING
  )
  
  op(x, NA_letter)
}
