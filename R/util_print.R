#' Print sq object
#' 
#' @description Prints input \code{\link{sq}} object in a human-friendly form.  
#' 
#' @details \code{Print} method is used by default in each case of calling the 
#' \code{\link{sq}} object with default parameters. 
#' Only by explicit calling the \code{print} method parameters can be changed. 
#'  
#' \code{Print} checks if the input \code{\link{sq}} object is cleaned and includes 
#' this information alongside with type in the printed message. On the right side of 
#' the sequence, in angle brackets, the length of each sequence is printed (e.q. "<9>").
#' 
#' If the \code{max_sequences} parameter is supplied, the desired number of sequences 
#' is printed and this information is included in a message (e.q. "printed 1 out of 3"). 
#' Only \code{max_sequences} value smaller than the number of sequences in object 
#' affects the function. The default value indicating how many sequences should 
#' be printed is 10, but it can be changed in \code{\link[=tidysq-options]{package options}}. 
#' 
#' Default value of \code{use_color} parameter is \code{TRUE} - sequences are printed
#' in green, while empty sequences, NA character and dots in gray. If this option is disabled, 
#' all sequences are in default color of console.
#' 
#' The \code{letters_sep} parameter indicates how the letters should be separated 
#' (they are not by default). Any character string can be supplied but 
#' \code{\link{NA_character_}}.
#' 
#' If sequences are too long, only leading characters are printed (as many as possible
#' in single line) and following dots indicating that sequence is truncated.
#' 
#' If sequences contain \code{\link{NA}} (‘Not Available’ / Missing Values) values, they 
#' are printed as "!" character, but it can be changed in 
#' \code{\link[=tidysq-options]{package options}}.
#' 
#' This is overloaded function from base package. It is selected when \code{\link{sq}} 
#' object is used as a parameter for print function. To see the generic function 
#' page, check \code{\link[base:print]{here}}.
#' 
#' @param x \code{\link{sq}} object
#' @param max_sequences \code{numeric} value indicating how many sequences 
#' should be printed
#' @param use_color \code{logical} value indicating if sequences should 
#' be colored
#' @param letters_sep \code{character} value indicating how the letters 
#' should be separated
#' @param ... 	further arguments passed to or from other methods. 
#' Unused.
#' 
#' @examples
#' 
#' # Creating sq objects using construct_sq:
#' sq_ami <- construct_sq(c("MIAANYTWIL","TIAALGNIIYRAIE", 
#'                          "NYERTGHLI", "MAYXXXIALN"), type = "ami")
#' sq_dna <- construct_sq(c("ATGCAGGA", "GACCGAACGAN", 
#'                          "TGACGAGCTTA"), type = "dna")
#' sq_unt <- construct_sq(c("ATGCAGGA!", "TGACGAGCTTA", "", "TIAALGNIIYRAIE"))
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
#' # Printing of the cleaned object:
#' clean(sq_dna)
#' print(clean(sq_dna), letters_sep = "-", use_color = FALSE)
#' 
#' @seealso \code{\link{sq}} \code{\link{clean}} \code{\link{tidysq-options}}
#' @name print.sq
NULL

#' @importFrom cli cli_text
#' @export
obj_print_header.sq <- function(x, ...) {
  if (length(.get_sq_type(x)) != 1)
    cli_text("sq (improper subtype!)")
  else
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
                                max_sequences = getOption("tidysq_p_max_sequences")) {
  if (length(x) > max_sequences)
    cli_text("printed {max_sequences} out of {length(x)}")
  invisible(x)
}

#' @importFrom cli col_green col_silver
#' @importFrom pillar pillar_shaft
#' @export
pillar_shaft.sq <- function(x, ...) {
  # color NA's
  na_character(alphabet(x)) <- col_silver(.get_na_char())
  
  .pillar_shaft_sq(x, "", col_green)
}

#' @importFrom cli col_cyan
#' @importFrom pillar pillar_shaft
#' @export
pillar_shaft.encsq <- function(x, ...) {
  alphabet(x) <- format(alphabet(x), digits = 1, scientific = FALSE)
  
  .pillar_shaft_sq(x, "\u00a0", col_cyan)
}
