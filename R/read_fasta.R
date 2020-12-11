#' Read a FASTA file
#'
#' @templateVar alph_null_ok TRUE
#'
#' @description Reads a FASTA file that contains nucleotide or amino acid
#' sequences and returns a \code{\link[tibble]{tibble}} with obtained data.
#'
#' @param file_name [\code{character(1)}]\cr
#'  Absolute path to file or url to read from.
#' @template alphabet
#' @template NA_letter
#' @template safe_mode
#' @template ignore_case
#'
#' @return A \code{\link[tibble]{tibble}} with number of rows equal to the
#' number of sequences and two columns:
#' \itemize{
#' \item{name}{specifies name of a sequence, used in functions like
#'  \code{\link{find_motifs}}}
#' \item{sq}{contains extracted sequence itself}
#' }
#'
#' @details
#' All rules of creating \code{sq} objects are the same as in \code{\link{sq}}.
#'
#' @examples
#' fasta_file <- system.file(package = "tidysq", "examples/example_aa.fasta")
#'
#' # In this case, these two calls are equivalent in result:
#' read_fasta(fasta_file)
#' read_fasta(fasta_file, alphabet = "ami_bsc")
#'
#' \dontrun{
#' # It's possible to read FASTA file from URL:
#' read_fasta("https://www.uniprot.org/uniprot/P28307.fasta")
#' }
#'
#' @family input_functions
#' @seealso \code{\link[base]{readLines}}
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
    #TODO: issue #56
    alphabet <- sq_alphabet(alphabet, "atp")
  }
  
  CPP_read_fasta(file_name, alphabet, NA_letter, ignore_case)
}