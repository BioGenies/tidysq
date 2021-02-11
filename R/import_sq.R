#' Import sq objects from other objects
#' 
#' @description Creates \code{\link[=sq-class]{sq}} object from object of class
#' from another package. Currently supported packages are \pkg{ape},
#' \pkg{bioseq}, \pkg{Bioconductor} and \pkg{seqinr}. For exact list of
#' supported classes and resulting types, see details.
#' 
#' @param object [\code{any(1)}]\cr
#'  An object of one of supported classes.
#' @template three-dots
#' 
#' @return A \code{\link[tibble]{tibble}} with \code{sq} column of
#' \code{\link[=sq-class]{sq}} type representing the same sequences as given
#' object; the object has a type corresponding to the input type; if given
#' sequences have names, output \code{\link[tibble]{tibble}} will also have
#' another column \code{name} with those names
#' 
#' @details
#' Currently supported classes are as follows:
#' \itemize{
#' \item \code{ape}:
#'  \itemize{
#'  \item \code{AAbin} - imported as \strong{ami_bsc}
#'  \item \code{DNAbin} - imported as \strong{dna_bsc}
#'  \item \code{alignment} - exact type is guessed within \code{\link{sq}}
#'   function
#'  }
#' \item \code{bioseq}:
#'  \itemize{
#'  \item \code{bioseq_aa} - imported as \strong{ami_ext}
#'  \item \code{bioseq_dna} - imported as \strong{dna_ext}
#'  \item \code{bioseq_rna} - imported as \strong{rna_ext}
#'  }
#' \item \code{Biostrings}:
#'  \itemize{
#'  \item \code{AAString} - imported as \strong{ami_ext} with exactly one
#'   sequence
#'  \item \code{AAStringSet} - imported as \strong{ami_ext}
#'  \item \code{DNAString} - imported as \strong{dna_ext} with exactly one
#'   sequence
#'  \item \code{DNAStringSet} - imported as \strong{dna_ext}
#'  \item \code{RNAString} - imported as \strong{rna_ext} with exactly one
#'   sequence
#'  \item \code{RNAStringSet} - imported as \strong{rna_ext}
#'  \item \code{BString} - imported as \strong{unt} with exactly one
#'   sequence
#'  \item \code{BStringSet} - imported as \strong{unt}
#'  \item \code{XStringSetList} - each element of a list can be imported as
#'   a separate \code{\link[tibble]{tibble}}, resulting in a list of tibbles;
#'   if passed argument \code{separate = FALSE}, these tibbles are bound into
#'   one bigger tibble
#'  }
#' \item \code{seqinr}:
#'  \itemize{
#'  \item \code{SeqFastaAA} - imported as \strong{ami_bsc}
#'  \item \code{SeqFastadna} - imported as \strong{dna_bsc}
#'  }
#' }
#'
#' Providing object of class other than specified will result in an error.
#'
#' @examples
#' # ape example
#' library(ape)
#' ape_dna <- as.DNAbin(list(one = c("C", "T", "C", "A"), two = c("T", "G", "A", "G", "G")))
#' import_sq(ape_dna)
#'
#' # bioseq example
#' library(bioseq)
#' bioseq_rna <- new_rna(c(one = "ANBRY", two = "YUTUGGN"))
#' import_sq(bioseq_rna)
#'
#' # Biostrings example
#' library(Biostrings)
#' Biostrings_ami <- AAStringSet(c(one = "FEAPQLIWY", two = "EGITENAK"))
#' import_sq(Biostrings_ami)
#'
#' # seqinr example
#' library(seqinr)
#' seqinr_dna <- as.SeqFastadna(c("C", "T", "C", "A"), name = "one")
#' import_sq(seqinr_dna)
#'
#' @family input_functions
#' @seealso \code{\link[=sq-class]{sq class}}
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
    }), alphabet = "ami_bsc", ...)
  } else if (is.list(x)) {
    x <- sq(vapply(x, function(i) {
      toupper(paste(i, collapse = ""))
    }, character(1)), alphabet = "ami_bsc", ...)
  } else if (is.character(x)) {
    # Sometimes obtained e.g. by extracting an element from whole AAbin list
    # Using code for list case should work as well, separate case is probably an overkill
    x <- sq(toupper(paste(x, collapse = "")), alphabet = "ami_bsc", ...)
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
    }), alphabet = "dna_bsc", ...)
  } else if (is.list(x)) {
    x <- sq(vapply(x, function(i) {
      toupper(paste(i, collapse = ""))
    }, character(1)), alphabet = "dna_bsc", ...)
  } else if (is.character(x)) {
    # Sometimes obtained e.g. by extracting an element from whole DNAbin list
    # Using code for list case should work as well, separate case is probably an overkill
    x <- sq(toupper(paste(x, collapse = "")), alphabet = "dna_bsc", ...)
  }
  bind_into_sqibble(x, labels(object))
}

#' @export
import_sq.alignment <- function(object, ...) {
  # From package `ape`
  bind_into_sqibble(sq(object[["seq"]], ...), object[["nam"]])
}

#' @export
import_sq.bioseq_aa <- function(object, ...) {
  # From package `bioseq`
  bind_into_sqibble(sq(as.character(object), alphabet = "ami_ext"), names(object))
}

#' @export
import_sq.bioseq_dna <- function(object, ...) {
  # From package `bioseq`
  bind_into_sqibble(sq(as.character(object), alphabet = "dna_ext"), names(object))
}

#' @export
import_sq.bioseq_rna <- function(object, ...) {
  # From package `bioseq`
  bind_into_sqibble(sq(as.character(object), alphabet = "rna_ext"), names(object))
}

#' @export
import_sq.AAString <- function(object, ...) {
  # From package `Biostrings`
  x <- sq(as.character(object), alphabet = "ami_ext", ...)
  bind_into_sqibble(x, names(object))
}

#' @export
import_sq.AAStringSet <- function(object, ...) {
  # From package `Biostrings`
  x <- sq(as.character(object), alphabet = "ami_ext", ...)
  bind_into_sqibble(x, names(object))
}

#' @export
import_sq.DNAString <- function(object, ...) {
  # From package `Biostrings`
  x <- sq(as.character(object), alphabet = "dna_ext", ...)
  bind_into_sqibble(x, names(object))
}

#' @export
import_sq.DNAStringSet <- function(object, ...) {
  # From package `Biostrings`
  x <- sq(as.character(object), alphabet = "dna_ext", ...)
  bind_into_sqibble(x, names(object))
}

#' @export
import_sq.RNAString <- function(object, ...) {
  # From package `Biostrings`
  x <- sq(as.character(object), alphabet = "rna_ext", ...)
  bind_into_sqibble(x, names(object))
}

#' @export
import_sq.RNAStringSet <- function(object, ...) {
  # From package `Biostrings`
  x <- sq(as.character(object), alphabet = "rna_ext", ...)
  bind_into_sqibble(x, names(object))
}

#' @export
import_sq.BString <- function(object, ...) {
  # From package `Biostrings`
  x <- sq(as.character(object), alphabet = "unt", ...)
  bind_into_sqibble(x, names(object))
}

#' @export
import_sq.BStringSet <- function(object, ...) {
  # From package `Biostrings`
  x <- sq(as.character(object), alphabet = "unt", ...)
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
  x <- sq(paste(object, collapse = ""), alphabet = "ami_bsc", ...)
  bind_into_sqibble(x, attr(object, "name"))
}

#' @export
import_sq.SeqFastadna <- function(object, ...) {
  # From package `seqinr`
  x <- sq(paste(object, collapse = ""), alphabet = "dna_bsc", ...)
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
