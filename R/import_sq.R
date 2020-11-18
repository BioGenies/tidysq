#' Import sq objects from other objects
#' 
#' Creates \code{\link[=sq]{sq object}} from \code{object} of class from another package.
#' Currently supported packages are \pkg{ape} with its formats (\code{AAbin} and \code{DNAbin}),
#' \pkg{Bioconductor} (\code{AAStringSet}, \code{DNAStringSet}) and
#' \pkg{seqinr} (\code{SeqFastaAA}, \code{SeqFastadna}).
#' 
#' @param object - an object of one of classes: \code{AAbin}, \code{DNAbin}, \code{AAStringSet}, 
#' \code{DNAStringSet}, \code{SeqFastaAA}, \code{SeqFastadna}.
#' @param ... - additional arguments passed to the function.
#' 
#' @return A \code{\link[tibble]{tibble}} with \code{sq} column of \code{\link{sq}} type 
#' representing the same 
#' sequences as given object; the object has a type corresponding to the input type; if given
#' sequences had names, output \code{\link[tibble]{tibble}} has also another column 
#' \code{name} with those names
#' 
#' @details 
#' Providing object of class other than specified will result in error.
#' 
#' @seealso \code{\link{export_sq}} \code{\link{sq}}
#' @export
import_sq <- function(object, ...)
  UseMethod("import_sq")

#' @export
import_sq.default <- function(object, ...)
  stop("import_sq() function cannot handle objects with this class", call. = FALSE)

#' @export
import_sq.AAbin <- function(object, ...) {
  # From package `ape`
  x <- as.character(object)
  if (is.matrix(x)) {
    x <- sq(apply(x, 1, function(i) {
      toupper(paste(i, collapse = ""))
    }), "ami_bsc")
  } else if (is.list(x)) {
    x <- sq(vapply(x, function(i) {
      toupper(paste(i, collapse = ""))
    }, character(1)), "ami_bsc")
  } else if (is.character(x)) {
    # Sometimes obtained e.g. by extracting an element from whole AAbin list
    # Using code for list case should work as well, separate case is probably an overkill
    x <- sq(toupper(paste(x, collapse = "")), "ami_bsc")
  }
  bind_into_sqibble(x, labels(object))
}

#' @export
import_sq.DNAbin <- function(object, ...) {
  # From package `ape`
  x <- as.character(object)
  if (is.matrix(x)) {
    x <- sq(apply(x, 1, function(i) {
      toupper(paste(i, collapse = ""))
    }), "dna_bsc")
  } else if (is.list(x)) {
    x <- sq(vapply(x, function(i) {
      toupper(paste(i, collapse = ""))
    }, character(1)), "dna_bsc")
  } else if (is.character(x)) {
    # Sometimes obtained e.g. by extracting an element from whole DNAbin list
    # Using code for list case should work as well, separate case is probably an overkill
    x <- sq(toupper(paste(x, collapse = "")), "dna_bsc")
  }
  bind_into_sqibble(x, labels(object))
}

#' @export
import_sq.alignment <- function(object, ...) {
  # From package `ape`
  bind_into_sqibble(sq(object[["seq"]]), object[["nam"]])
}

#' @export
import_sq.AAString <- function(object, ...) {
  # From package `Biostrings`
  x <- sq(as.character(object), alphabet = "ami_ext")
  bind_into_sqibble(x, names(object))
}

#' @export
import_sq.AAStringSet <- function(object, ...) {
  # From package `Biostrings`
  x <- sq(as.character(object), alphabet = "ami_ext")
  bind_into_sqibble(x, names(object))
}

#' @export
import_sq.DNAString <- function(object, ...) {
  # From package `Biostrings`
  x <- sq(as.character(object), alphabet = "dna_ext")
  bind_into_sqibble(x, names(object))
}

#' @export
import_sq.DNAStringSet <- function(object, ...) {
  # From package `Biostrings`
  x <- sq(as.character(object), alphabet = "dna_ext")
  bind_into_sqibble(x, names(object))
}

#' @export
import_sq.RNAString <- function(object, ...) {
  # From package `Biostrings`
  x <- sq(as.character(object), alphabet = "rna_ext")
  bind_into_sqibble(x, names(object))
}

#' @export
import_sq.RNAStringSet <- function(object, ...) {
  # From package `Biostrings`
  x <- sq(as.character(object), alphabet = "rna_ext")
  bind_into_sqibble(x, names(object))
}

#' @export
import_sq.BString <- function(object, ...) {
  # From package `Biostrings`
  x <- sq(as.character(object), alphabet = "unt")
  bind_into_sqibble(x, names(object))
}

#' @export
import_sq.BStringSet <- function(object, ...) {
  # From package `Biostrings`
  x <- sq(as.character(object), alphabet = "unt")
  bind_into_sqibble(x, names(object))
}

#' @export
import_sq.XStringSetList <- function(object, ...) {
  # From package `Biostrings`
  # Allowing the user to choose whether import should preserve there being separate sets
  # Is not a list even though behaves like a list (for us), so can't rely on this dispatch
  import_sq.list(object, ...)
}

#' @export
import_sq.SeqFastaAA <- function(object, ...) {
  # From package `seqinr`
  x <- sq(paste(object, collapse = ""), "ami_bsc")
  bind_into_sqibble(x, attr(object, "name"))
}

#' @export
import_sq.SeqFastadna <- function(object, ...) {
  # From package `seqinr`
  x <- sq(paste(object, collapse = ""), "dna_bsc")
  bind_into_sqibble(x, attr(object, "name"))
}

#' @importFrom dplyr bind_rows
#' @export
import_sq.list <- function(object, separate = TRUE, ...) {
  # Whenever a class is attributed to each element and not to the whole list
  # For example - `seqinr` and its SeqFastadna
  if (separate)
    lapply(object, import_sq)
  else
    do.call(bind_rows, lapply(object, import_sq))
}

bind_into_sqibble <- function(x, name = NULL) {
  if (is.null(name))
    tibble(sq = x)
  else 
    tibble(name = name, sq = x)
}
