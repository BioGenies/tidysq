#' @param safe_mode [\code{logical(1)}]\cr
#'  Default value is \code{FALSE}. When turned on, safe mode guarantees that
#'  \code{NA} appears within a sequence if and only if input sequence contains
#'  value passed with \code{NA_letter}. This means that resulting type might be
#'  different to the one passed as argument, if there are letters in a sequence
#'  that does not appear in the original alphabet.
