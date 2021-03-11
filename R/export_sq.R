#' Export sq objects into other formats
#'
#' @templateVar name_null_ok TRUE
#' 
#' @description Converts object of class \code{\link[=sq-class]{sq}} to a class
#' from another package. Currently supported packages are \pkg{ape},
#' \pkg{bioseq}, \pkg{Bioconductor} and \pkg{seqinr}. For exact list of
#' supported classes and resulting types, see details.
#'
#' @template x
#' @param export_format [\code{character(1)}]\cr
#'  A string indicating desired class (with specified package for clarity).
#' @template name
#' @template three-dots
#'
#' @return An object with the format specified in the parameter. To find 
#' information about the detailed structure of this object, see documentation 
#' of these objects.
#'
#' @details
#' Currently supported formats are as follows (grouped by \code{sq} types):
#' \itemize{
#' \item \strong{ami}:
#'  \itemize{
#'  \item \code{"ape::AAbin"}
#'  \item \code{"bioseq::bioseq_aa"}
#'  \item \code{"Biostrings::AAString"}
#'  \item \code{"Biostrings::AAStringSet"}
#'  \item \code{"seqinr::SeqFastaAA"}
#'  }
#' \item \strong{dna}:
#'  \itemize{
#'  \item \code{"ape::DNAbin"}
#'  \item \code{"bioseq::bioseq_dna"}
#'  \item \code{"Biostrings::DNAString"}
#'  \item \code{"Biostrings::DNAStringSet"}
#'  \item \code{"seqinr::SeqFastadna"}
#'  }
#' \item \strong{rna}:
#'  \itemize{
#'  \item \code{"bioseq::bioseq_rna"}
#'  \item \code{"Biostrings::RNAString"}
#'  \item \code{"Biostrings::RNAStringSet"}
#'  }
#' }
#'
#' @examples
#' # DNA and amino acid sequences can be exported to most packages
#' sq_ami <- sq(c("MVVGL", "LAVPP"), alphabet = "ami_bsc")
#' export_sq(sq_ami, "ape::AAbin")
#' export_sq(sq_ami, "bioseq::bioseq_aa")
#' export_sq(sq_ami, "Biostrings::AAStringSet", c("one", "two"))
#' export_sq(sq_ami, "seqinr::SeqFastaAA")
#'
#' sq_dna <- sq(c("TGATGAAGCGCA", "TTGATGGGAA"), alphabet = "dna_bsc")
#' export_sq(sq_dna, "ape::DNAbin", name = c("one", "two"))
#' export_sq(sq_dna, "bioseq::bioseq_dna")
#' export_sq(sq_dna, "Biostrings::DNAStringSet")
#' export_sq(sq_dna, "seqinr::SeqFastadna")
#'
#' # RNA sequences are limited to Biostrings and bioseq
#' sq_rna <- sq(c("NUARYGCB", "", "DRKCNYBAU"), alphabet = "rna_ext")
#' export_sq(sq_rna, "bioseq::bioseq_rna")
#' export_sq(sq_rna, "Biostrings::RNAStringSet")
#'
#' # Biostrings can export single sequences to simple strings as well
#' export_sq(sq_dna[1], "Biostrings::DNAString")
#'
#' @family output_functions
#' @seealso \code{\link[=sq-class]{sq class}}
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
    `bioseq::bioseq_aa` = {
      assert_package_installed("bioseq")
      bioseq::new_aa(setNames(as.character(x), name))
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
          seqinr::as.SeqFastaAA(sequence, name = seq_name)
        }, unpack(x, "STRINGS"), name, SIMPLIFY = FALSE, USE.NAMES = FALSE)
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
    `bioseq::bioseq_dna` = {
      assert_package_installed("bioseq")
      bioseq::new_dna(setNames(as.character(x), name))
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
          seqinr::as.SeqFastadna(sequence, name = seq_name)
        }, unpack(x, "STRINGS"), name, SIMPLIFY = FALSE, USE.NAMES = FALSE)
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
    `bioseq::bioseq_rna` = {
      assert_package_installed("bioseq")
      bioseq::new_rna(setNames(as.character(x), name))
    },
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
