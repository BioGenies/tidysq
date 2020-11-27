#' Remove sequences that contain NA values
#' 
#' Removes sequences containing ambiguous elements or removes \code{\link[=sq-class]{NA values}}
#' from sequences in a \code{\link[=sq-class]{sq}} object.
#' 
#' @inheritParams remove_ambiguous
#' @template NA_letter
#' @template three-dots
#'  
#' @return A \code{\link[=sq-class]{sq}} object with the same type as input type. Sequences not containing
#' any \code{\link[=sq-class]{NA}} values are left untouched.
#' 
#' @details This function allows removal of sequences containing \code{\link[=sq-class]{NA}} values.
#' By default, whole sequences containing ambiguous elements are removed 
#' and \code{\link[=sq-class]{NULL}} (empty) sequences are introduced in their place. If
#' \code{only_elements = TRUE} then only \code{\link[=sq-class]{NA}} values are removed
#' from sequences in \code{sq} object. \code{\link[=sq-class]{NULL}} values (empty sequences)
#' can be identified using \code{\link{is_empty_sq}} function.
#' 
#' \code{NA} may be introduced as a result of using functions like 
#' \code{\link{substitute_letters}} or \code{\link{bite}}. They also appear in sequences if
#' you are reading file using \code{\link{read_fasta}} or constructing \code{sq} object from
#' \code{\link{character}} vector with \code{\link{construct_sq}} in 
#' \code{\link[=fast-mode]{fast mode}} and there are letters in file or in strings other than
#' specified.
#'
#' @seealso \code{\link[=sq-class]{sq}} \code{\link{is_empty_sq}} \code{\link{substitute_letters}}
#' \code{\link{bite}}
#' @export
remove_na <- function(x, by_letter = FALSE, ...) {
  assert_flag(by_letter)
  
  UseMethod("remove_na")
}

#' @export
remove_na.default <- function(x, by_letter = FALSE, ...)
  stop("'remove_na' isn't implemented for this type of object", call. = FALSE)

#' @rdname remove_na
#' @export
remove_na.sq <- function(x, by_letter = FALSE, ...,
                         NA_letter = getOption("tidysq_NA_letter")) {
  assert_string(NA_letter, min.chars = 1)
  
  CPP_remove_NA(x, by_letter, NA_letter)
}
