#' @export
read_fasta <- function(file_name,
                       alphabet = NULL,
                       NA_letter = getOption("tidysq_NA_letter"),
                       safe_mode = getOption("tidysq_safe_mode"),
                       ignore_case = FALSE) {
  assert_character(file_name, any.missing = FALSE)
  assert_flag(safe_mode)
  assert_string(NA_letter)
  assert_character(alphabet, any.missing = FALSE, min.len = 0, unique = TRUE, null.ok = TRUE)
  assert_flag(ignore_case)
  
  if (is.null(alphabet)) {
    alphabet <- CPP_sample_fasta(file_name, if (safe_mode) Inf else 4096, 
                                 NA_letter, ignore_case)
    alphabet <- guess_standard_alphabet(alphabet)
  } else if (length(alphabet) == 1) {
    type <- interpret_type(alphabet)
    if (type == "unt") {
      alphabet <- CPP_sample_fasta(file_name, Inf, NA_letter, ignore_case)
    } else {
      alphabet <- get_standard_alphabet(type)
      if (safe_mode) {
        actual_alphabet <- CPP_sample_fasta(file_name, Inf, NA_letter, ignore_case)
        if (!identical(actual_alphabet, alphabet)){
          warning("Detected letters that do not match specified type!")
          alphabet <- actual_alphabet
        }
      }
    }
  } else {
    #TODO: safe mode should also be implemented for atp
    alphabet <- sq_alphabet(alphabet, "atp")
  }
  
  CPP_read_fasta(file_name, alphabet, NA_letter, ignore_case)
}