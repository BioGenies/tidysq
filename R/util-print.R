# TODO: issue #59

#' @importFrom cli cli_text
#' @export
obj_print_header.sq <- function(x, ...) {
  cli_text("{vec_ptype_full(x)} sequences list",
           if (length(x) == 0) " of length 0" else ":")
}

#' @importFrom cli cat_line
#' @export
obj_print_data.sq <- function(x, ...) {
  # for some unknown reason had to copy this code from original obj_print_data
  if (length(x) == 0) {
    return()
  }
  cat_line(format(x, ...))
  invisible(x)
}

#' @importFrom cli cli_text
#' @export
obj_print_footer.sq <- function(x, ...,
                                max_sequences = getOption("tidysq_print_max_sequences")) {
  if (length(x) > max_sequences)
    cli_text("printed {max_sequences} out of {length(x)}")
  invisible(x)
}

#' @importFrom cli col_green
#' @importFrom pillar pillar_shaft
#' @export
pillar_shaft.sq <- function(x, ...,
                            NA_letter = getOption("tidysq_NA_letter"),
                            max_pillar_width = getOption("tidysq_pillar_max_width")) {
  assert_string(NA_letter)
  assert_integerish(max_pillar_width, lower = 6, len = 1)

  pillar_shaft_sq(x, "", NA_letter, col_green, max_pillar_width)
}

#' @importFrom cli col_cyan
#' @importFrom pillar pillar_shaft
#' @export
pillar_shaft.encsq <- function(x, ...,
                               NA_letter = getOption("tidysq_NA_letter"),
                               max_pillar_width = getOption("tidysq_pillar_max_width")) {
  assert_string(NA_letter)
  assert_integerish(max_pillar_width, lower = 6, len = 1)

  alphabet(x) <- format(alphabet(x), digits = 1, scientific = FALSE)
  
  pillar_shaft_sq(x, "\u00a0", NA_letter, col_cyan, max_pillar_width)
}

#' Print sq object
#' 
#' @description Prints input \code{\link[=sq-class]{sq}} object in a
#' human-friendly form.
#'
#' @template x
#' @param max_sequences [\code{integer(1)}]\cr
#'  How many sequences should be printed.
#' @param use_color [\code{logical(1)}]\cr
#'  Should sequences be colored?
#' @param letters_sep [\code{character(1)}]\cr
#'  How the letters should be separated.
#' @template NA_letter
#' @template three-dots
#' 
#' @return An object that was passed as the first argument to the function.
#' It is returned invisibly (equivalent of \code{invisible(x)})
#' 
#' @details
#' \code{print} method is often called implicitly by calling variable name.
#' Only explicit calling of this method allows its parameters to be changed.
#'
#' Printed information consists of three parts:
#' \itemize{
#' \item First line is always a header that contains info about the type of
#'  sequences contained.
#' \item The next part is the content. Each sequence has its own line, but not
#'  all sequences are printed. The number of printed sequences is limited by
#'  parameter \code{max_sequences}, defaulting to 10. These sequences are
#'  printed with:
#'  \itemize{
#'  \item left-aligned index of sequence in square brackets (e.g. \code{[3]}),
#'  \item left-aligned sequence data (more about it in paragraph below),
#'  \item right-aligned sequence length in angle brackets (e.g. \code{<27>}).
#'  }
#' \item Finally, if number of sequences is greater than \code{max_sequences},
#'  then a footer is displayed with how many sequences are there and how many
#'  were printed.
#' }
#'
#' Each sequence data is printed as letters. If sequence is too long to fit in
#' one line, then only a subsequence is displayed - a subsequence that begins
#' from the first letter. Sequence printing is controlled by \code{letters_sep}
#' and \code{NA_letter} parameters. The first one specifies a string that should
#' be inserted between any two letters. By default it's empty when all letters
#' are one character in length; and a space otherwise. \code{NA_letter} dictates
#' how \code{NA} values are displayed, by default it's an exclamation mark
#' ("\code{!}").
#'
#' Most consoles support color printing, but when any of these do not, then the
#' user might use \code{use_color} parameter set to \code{FALSE} - or better
#' yet, change related option value, where said option is called
#' \code{"tidysq_print_use_color"}.
#'
#' @examples
#' # Creating objects to work on:
#' sq_ami <- sq(c("MIAANYTWIL","TIAALGNIIYRAIE", "NYERTGHLI", "MAYXXXIALN"),
#'              alphabet = "ami_ext")
#' sq_dna <- sq(c("ATGCAGGA", "GACCGNBAACGAN", "TGACGAGCTTA"),
#'              alphabet = "dna_bsc")
#' sq_unt <- sq(c("ATGCAGGA?", "TGACGAGCTTA", "", "TIAALGNIIYRAIE"))
#'
#' # Printing without explicit function calling with default parameters:
#' sq_ami
#' sq_dna
#' sq_unt
#'
#' # Printing with explicit function calling and specific parameters:
#' print(sq_ami)
#' print(sq_dna, max_sequences = 1, use_color = FALSE)
#' print(sq_unt, letters_sep = ":")
#'
#' @family display_functions
#' @name sqprint
#' @aliases sq-print
NULL
