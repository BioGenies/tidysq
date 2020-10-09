#' Export sq objects into other formats
#' 
#' Convert object of class \code{\link{sq}} to another class from another package. Currently 
#' supported packages are \pkg{ape} with its formats (\code{AAbin} and \code{DNAbin}),
#' \pkg{Bioconductor} (\code{AAStringSet}, \code{DNAStringSet}) and
#' \pkg{seqinr} (\code{SeqFastaAA}, \code{SeqFastadna}).
#' @inheritParams write_fasta
#' @param export_format a \code{\link{character}} string indicating package and the destination 
#' class; it should be one of the following: "seqinr::SeqFastaAA", "ape::AAbin", 
#' "Biostrings::AAStringSet", "seqinr::SeqFastadna", "ape::DNAbin", "Biostrings::DNAStringSet".
#' @param ... - additional arguments passed to the function.
#' 
#' @examples 
#' sq_ami <- construct_sq(c("MVVGL", "LAVPP"))
#' export_sq(sq_ami, "ape::AAbin")
#' export_sq(sq_ami, "Biostrings::AAStringSet", c("one", "two"))
#' export_sq(sq_ami, "seqinr::SeqFastaAA")
#' 
#' sq_dna <- construct_sq(c("TGATGAAGCGCA", "TTGATGGGAA"))
#' export_sq(sq_dna, "ape::DNAbin", name = c("one", "two"))
#' export_sq(sq_dna, "Biostrings::DNAStringSet")
#' export_sq(sq_dna, "seqinr::SeqFastadna")
#' @seealso \code{\link{sq}} \code{\link{import_sq}}
#' @export
export_sq <- function(x, export_format, name = NULL, ...) {
  assert_string(export_format)
  assert_character(name, len = vec_size(x), null.ok = TRUE)
  
  UseMethod("export_sq")
}

#' @export
export_sq.default <- function(x, export_format, name = NULL, ...)
  stop("export_sq() function cannot export objects of this class", call. = FALSE)

#' @export
export_sq.amisq <- function(x, export_format, name = NULL, ...) {
  switch (export_format,
    `ape::AAbin` = {
      .check_is_installed("ape")
      ape::as.AAbin(setNames(lapply(.unpack_from_sq(x, "char"), `attributes<-`, NULL), name))
    },
    `Biostrings::AAString` = {
      .check_is_installed("Biostrings")
      if (vec_size(x) != 1)
        stop("sq object must contain exactly one sentence; otherwise use \"Biostrings::AAStringSet\"", call. = FALSE)
      Biostrings::AAString(setNames(unlist(.unpack_from_sq(x, "string")), name))
    },
    `Biostrings::AAStringSet` = {
      .check_is_installed("Biostrings")
      Biostrings::AAStringSet(setNames(unlist(.unpack_from_sq(x, "string")), name))
    },
    `seqinr::SeqFastaAA` = {
      .check_is_installed("seqinr")
      if (is.null(name)) {
        lapply(.unpack_from_sq(x, "char"), seqinr::as.SeqFastaAA)
      } else {
        mapply(function(sequence, seq_name) {
          `attr<-`(seqinr::as.SeqFastaAA(sequence), "name", seq_name)
        }, .unpack_from_sq(x, "char"), name, SIMPLIFY = FALSE)
      }
    },
    {
      stop("exporting to this format is not yet supported; else, maybe you misspelled export_format parameter?", call. = FALSE)
    }
  )
  
}

#' @export
export_sq.dnasq <- function(x, export_format, name = NULL, ...) {
  switch (export_format,
    `ape::DNAbin` = {
      .check_is_installed("ape")
      ape::as.DNAbin(setNames(lapply(.unpack_from_sq(x, "char"), `attributes<-`, NULL), name))
    },
    `Biostrings::DNAString` = {
      .check_is_installed("Biostrings")
      if (vec_size(x) != 1)
        stop("sq object must contain exactly one sentence; otherwise use \"Biostrings::DNAStringSet\"", call. = FALSE)
      Biostrings::DNAString(setNames(unlist(.unpack_from_sq(x, "string")), name))
    },
    `Biostrings::DNAStringSet` = {
      .check_is_installed("Biostrings")
      Biostrings::DNAStringSet(setNames(unlist(.unpack_from_sq(x, "string")), name))
    },
    `seqinr::SeqFastadna` = {
      .check_is_installed("seqinr")
      if (is.null(name)) {
        lapply(.unpack_from_sq(x, "char"), seqinr::as.SeqFastadna)
      } else {
        mapply(function(sequence, seq_name) {
          `attr<-`(seqinr::as.SeqFastadna(sequence), "name", seq_name)
        }, .unpack_from_sq(x, "char"), name, SIMPLIFY = FALSE)
      }
    },
    {
      stop("exporting to this format is not yet supported; else, maybe you misspelled export_format parameter?", call. = FALSE)
    }
  )
}

#' @export
export_sq.rnasq <- function(x, export_format, name = NULL, ...) {
  switch (export_format,
    `Biostrings::RNAString` = {
      .check_is_installed("Biostrings")
      if (vec_size(x) != 1)
        stop("sq object must contain exactly one sentence; otherwise use \"Biostrings::RNAStringSet\"", call. = FALSE)
      Biostrings::RNAString(setNames(unlist(.unpack_from_sq(x, "string")), name))
    },
    `Biostrings::RNAStringSet` = {
      .check_is_installed("Biostrings")
      Biostrings::RNAStringSet(setNames(unlist(.unpack_from_sq(x, "string")), name))
    },
    {
      stop("exporting to this format is not yet supported; else, maybe you misspelled export_format parameter?", call. = FALSE)
    }
  )
}
