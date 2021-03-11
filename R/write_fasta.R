#' Save sq to fasta file
#'
#' @templateVar name_null_ok FALSE
#' 
#' @description Writes \code{\link[=sq-class]{sq}} objects with their names to
#' a fasta file.
#'
#' @template x
#' @template name
#' @param file [\code{character(1)}]\cr
#'  Absolute path to file to write to.
#' @param width [\code{integer(1)}]\cr
#'  Maximum number of characters to put in each line of file. Must be positive.
#' @template NA_letter
#'
#' @return No value is returned.
#' 
#' @details
#' Whenever a name has more letters than \code{width} parameter, nothing
#' happens, as only sequences are split to fit within designated space.
#'
#' @examples
#' \dontrun{
#' sq_dna <- sq(c("ACTGCTG", "CTTAGA", "CCCT", "CTGAATGT"),
#'              alphabet = "dna_bsc")
#' write_fasta(sq_dna,
#'             c("bat", "cat", "rat", "elephant_swallowed_by_A_snake"),
#'             "~/fasta_rubbish/example.fasta")
#' }
#'
#' @family output_functions
#' @export
write_fasta <- function(x, name, file,
                        width = 80,
                        NA_letter = getOption("tidysq_NA_letter")) {
  assert_class(x, "sq")
  assert_character(name, len = vec_size(x), any.missing = FALSE)
  assert_string(file)
  assert_count(width, positive = TRUE)
  assert_string(NA_letter, min.chars = 1)
  
  CPP_write_fasta(x, name, file, width, NA_letter)
}