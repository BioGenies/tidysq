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
#' @template on_warning
#' @template ignore_case
#'
#' @return A \code{\link[tibble]{tibble}} with number of rows equal to the
#' number of sequences and two columns:
#' \itemize{
#' \item `name` -- specifies name of a sequence, used in functions like \code{\link{find_motifs}}
#' \item `sq` -- specifies name of a sequence, used in functions like \code{\link{find_motifs}}
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
                       on_warning = getOption("tidysq_on_warning"),
                       ignore_case = FALSE) {
  sq_from_source(file_name, alphabet, CPP_sample_fasta, CPP_read_fasta,
                 NA_letter, safe_mode, on_warning, ignore_case)
}
