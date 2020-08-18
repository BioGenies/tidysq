#' @importFrom pillar type_sum
#' @export
type_sum.amisq <- function(x) {
  "ami"
}

#' @importFrom cli cat_line
#' @export
type_sum.dnasq <- function(x) {
  "dna"
}

#' @importFrom pillar type_sum
#' @export
type_sum.rnasq <- function(x) {
  "rna"
}

#' @importFrom pillar type_sum
#' @export
type_sum.untsq <- function(x) {
  "unt"
}

#' @importFrom pillar type_sum
#' @export
type_sum.atpsq <- function(x) {
  "atp"
}

#' @importFrom pillar type_sum
#' @export
type_sum.encsq <- function(x) {
  "enc"
}

#' @importFrom pillar type_sum
#' @export
type_sum.clnsq <- function(x) {
  paste0("(c)", NextMethod())
}

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
#' 
#' @importFrom crayon col_nchar
#' @importFrom crayon green
#' @importFrom crayon silver
#' @export
print.sq <- function(x,  
                     max_sequences = getOption("tidysq_p_max_sequences"),
                     use_color = getOption("tidysq_p_use_color"),
                     letters_sep = NULL, ...) {
  .check_integer(max_sequences, "'max_sequences'")
  .check_logical(use_color, "'use_color'")
  .check_character(letters_sep, "'letters_sep'", single_elem = TRUE, 
                   allow_zero_len = TRUE, allow_null = TRUE)
  
  # color NA's
  alph <- .get_alph(x)
  na_char <- if (.get_color_opt()) silver(.get_na_char()) else .get_na_char()
  na_val <- .get_na_val(alph)
  alph[na_val] <- na_char
  
  # if parameter is NULL and all letters are length one, no space
  if (is.null(letters_sep)) {
    letters_sep <- if (all(col_nchar(alph[!is.na(alph)]) == 1)) "" else " "
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
#' @name sqprint
#' @aliases sq-print
NULL
