#' Save sq to fasta file
#'
#' @templateVar name_null_ok FALSE
#' 
#' Writes \code{\link[=sq-class]{sq}} objects with their names to a fasta file.
#'
#' @template x
#' @template name
#' @param file [\code{character(1)}]\cr
#'  Absolute path to file to write to.
#' @param width [\code{integer(1)}]\cr
#'  Maximum number of characters to put in each line of file. Must be positive.
#' @template NA_letter
#'
#' @export
write_fasta <- function(x, names, file, 
                        width = 80,
                        NA_letter = getOption("tidysq_NA_letter")) {
  assert_class(x, "sq")
  assert_character(names, len = vec_size(x), any.missing = FALSE)
  assert_string(file)
  assert_count(width, positive = TRUE)
  assert_string(NA_letter, min.chars = 1)
  
  CPP_write_fasta(x, names, file, width, NA_letter)
}