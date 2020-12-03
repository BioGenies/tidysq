#' Subset sequences from sq objects
#' 
#' @description Extracts a defined range of elements (amino acids or nucleotides) 
#' from a sequence.
#' 
#' @template x
#' @param indices [\code{integer}]\cr
#'  Indices to extract from each sequence. The function follows the normal R
#'  conventions for indexing vectors, including negative indices.
#' @template NA_letter
#' @template on_warning
#' 
#' @return \code{\link[=sq-class]{sq}} object of the same type as input sq, where each
#' element is a subsequence created by indexing corresponding sequence from 
#' input sq object with input indices.
#' 
#' @details
#' Amino acids and nucleic acid sequences are represented as \code{\link[=sq-class]{sq}}
#' object in the \code{\link{tidysq}} package. Often one needs to get only a 
#' single letter, or the sequence of a defined range from the original sequences. 
#' A subsequence is a sequence that can be derived from the original sequence 
#' by trimming some elements (letters) without changing the order of the 
#' remaining elements. To get a subsequence from each sequence contained in 
#' the \code{\link[=sq-class]{sq}} object with the same indices. This is for example
#' useful to extract a user-defined region from a sequence. 
#' 
#' The usage of \code{bite} follows the normal R conventions. For details 
#' refer to the R documentation (see 
#' \url{https://cran.r-project.org/doc/manuals/r-release/R-intro.html#Index-vectors}). 
#' Extracting indices not present in the sequence results in introducing 
#' \code{\link{NA}} (‘Not Available’ / Missing Value) values. 
#' Information about it is printed on a console depending on the value of option 
#' 'tidysq_a_bite_na' - it can be either a warning (default), an error,
#' a message or no information (you can check details in \code{\link{tidysq-options})}. 
#' \code{NA} values can be removed by using \code{\link{remove_na}} function.
#'
#' @family order_functions
#' @seealso \code{\link{remove_na}}
#' @export
bite <- function(x, indices, ...)
  UseMethod("bite")

#' @export
bite.default <- function(x, indices, ...)
  stop("method 'bite()' isn't implemented for this type of object", call. = FALSE)

#' @rdname bite
#' @export
bite.sq <- function(x, indices, ...,
                    NA_letter = getOption("tidysq_NA_letter"),
                    on_warning = getOption("tidysq_on_warning")) {
  assert_string(NA_letter, min.chars = 1)
  assert_warning_handling(on_warning)
  assert_integerish(indices, any.missing = FALSE, null.ok = TRUE)
  
  ret <- CPP_bite(x, indices, NA_letter)
  handle_warning_message(ret[["warning"]], on_warning)
  ret[["sq"]]
}
