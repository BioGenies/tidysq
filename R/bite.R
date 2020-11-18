#' Subset sequences from sq objects
#' 
#' @description Extracts a defined range of elements (amino acids or nucleotides) 
#' from a sequence.
#' 
#' @inheritParams reverse
#' @param indices a \code{\link{numeric}} vector of subsequence indices to extract from
#' each sequence. The function follows the normal R conventions for indexing 
#' vectors, including negative indices.
#' 
#' @return \code{\link{sq}} object of the same type as input sq, where each 
#' element is a subsequence created by indexing corresponding sequence from 
#' input sq object with input indices.
#' 
#' @details Amino acids and nucleic acid sequences are represented as \code{\link{sq}} 
#' object in the \code{\link{tidysq}} package. Often one needs to get only a 
#' single letter, or the sequence of a defined range from the original sequences. 
#' A subsequence is a sequence that can be derived from the original sequence 
#' by trimming some elements (letters) without changing the order of the 
#' remaining elements. To get a subsequence from each sequence contained in 
#' the \code{\link{sq}} object with the same indices. This is for example 
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
#' @seealso \code{\link{sq}} \code{\link{remove_na}} \code{\link{tidysq-options}}
#' @export
bite <- function(x, indices, ...,
                 NA_letter = getOption("tidysq_NA_letter"))
  UseMethod("bite")

#' @export
bite.default <- function(x, indices, ...,
                         NA_letter = getOption("tidysq_NA_letter"))
  stop("method 'bite()' isn't implemented for this type of object", call. = FALSE)

#' @export
bite.sq <- function(x, indices, ...,
                    NA_letter = getOption("tidysq_NA_letter")) {
  assert_string(NA_letter, min.chars = 1)
  assert_integerish(indices, any.missing = FALSE, null.ok = TRUE)
  
  ret <- CPP_bite(x, indices, NA_letter)
  if (ret[["warning"]] != "")
    .handle_opt_txt("tidysq_a_bite_na", ret[["warning"]])
  ret[["sq"]]
}
