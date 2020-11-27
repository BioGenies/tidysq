#' Export sq objects into other formats
#' 
#' Convert object of class \code{\link[=sq-class]{sq}} to another class from another package. Currently
#' supported packages are \pkg{ape} with its formats (\code{AAbin} and \code{DNAbin}),
#' \pkg{Bioconductor} (\code{AAStringSet}, \code{DNAStringSet}) and
#' \pkg{seqinr} (\code{SeqFastaAA}, \code{SeqFastadna}).
#' @inheritParams write_fasta
#' @param export_format a \code{\link{character}} string indicating package and the destination 
#' class; it should be one of the following: "seqinr::SeqFastaAA", "ape::AAbin", 
#' "Biostrings::AAStringSet", "seqinr::SeqFastadna", "ape::DNAbin", "Biostrings::DNAStringSet".
#' @template three-dots
#' 
#' @seealso \code{\link[=sq-class]{sq}} \code{\link{import_sq}}
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
export_sq.sq_ami_bsc <- function(x, export_format, name = NULL, ...) {
  switch (export_format,
    `ape::AAbin` = {
      assert_package_installed("ape")
      ape::as.AAbin(setNames(lapply(unpack(x, "STRINGS"), `attributes<-`, NULL), name))
    },
    `Biostrings::AAString` = {
      assert_package_installed("Biostrings")
      if (vec_size(x) != 1)
        stop("sq object must contain exactly one sentence; otherwise use \"Biostrings::AAStringSet\"", call. = FALSE)
      Biostrings::AAString(setNames(unlist(unpack(x, "STRING")), name))
    },
    `Biostrings::AAStringSet` = {
      assert_package_installed("Biostrings")
      Biostrings::AAStringSet(setNames(unlist(unpack(x, "STRING")), name))
    },
    `seqinr::SeqFastaAA` = {
      assert_package_installed("seqinr")
      if (is.null(name)) {
        lapply(unpack(x, "STRINGS"), seqinr::as.SeqFastaAA)
      } else {
        mapply(function(sequence, seq_name) {
          `attr<-`(seqinr::as.SeqFastaAA(sequence), "name", seq_name)
        }, unpack(x, "STRINGS"), name, SIMPLIFY = FALSE)
      }
    },
    {
      stop("exporting to this format is not yet supported; else, maybe you misspelled export_format parameter?", call. = FALSE)
    }
  )
}

#' @export
export_sq.sq_ami_ext <- export_sq.sq_ami_bsc

#' @export
export_sq.sq_dna_bsc <- function(x, export_format, name = NULL, ...) {
  switch (export_format,
    `ape::DNAbin` = {
      assert_package_installed("ape")
      ape::as.DNAbin(setNames(lapply(unpack(x, "STRINGS"), `attributes<-`, NULL), name))
    },
    `Biostrings::DNAString` = {
      assert_package_installed("Biostrings")
      if (vec_size(x) != 1)
        stop("sq object must contain exactly one sentence; otherwise use \"Biostrings::DNAStringSet\"", call. = FALSE)
      Biostrings::DNAString(setNames(unlist(unpack(x, "STRING")), name))
    },
    `Biostrings::DNAStringSet` = {
      assert_package_installed("Biostrings")
      Biostrings::DNAStringSet(setNames(unlist(unpack(x, "STRING")), name))
    },
    `seqinr::SeqFastadna` = {
      assert_package_installed("seqinr")
      if (is.null(name)) {
        lapply(unpack(x, "STRINGS"), seqinr::as.SeqFastadna)
      } else {
        mapply(function(sequence, seq_name) {
          `attr<-`(seqinr::as.SeqFastadna(sequence), "name", seq_name)
        }, unpack(x, "STRINGS"), name, SIMPLIFY = FALSE)
      }
    },
    {
      stop("exporting to this format is not yet supported; else, maybe you misspelled export_format parameter?", call. = FALSE)
    }
  )
}

#' @export
export_sq.sq_dna_ext <- export_sq.sq_dna_bsc

#' @export
export_sq.sq_rna_bsc <- function(x, export_format, name = NULL, ...) {
  switch (export_format,
    `Biostrings::RNAString` = {
      assert_package_installed("Biostrings")
      if (vec_size(x) != 1)
        stop("sq object must contain exactly one sentence; otherwise use \"Biostrings::RNAStringSet\"", call. = FALSE)
      Biostrings::RNAString(setNames(unlist(unpack(x, "STRING")), name))
    },
    `Biostrings::RNAStringSet` = {
      assert_package_installed("Biostrings")
      Biostrings::RNAStringSet(setNames(unlist(unpack(x, "STRING")), name))
    },
    {
      stop("exporting to this format is not yet supported; else, maybe you misspelled export_format parameter?", call. = FALSE)
    }
  )
}

#' @export
export_sq.sq_rna_ext <- export_sq.sq_rna_bsc
