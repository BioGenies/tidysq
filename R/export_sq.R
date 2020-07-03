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
#' 
#' @examples 
#' sq_ami <- construct_sq(c("MVVGL", "LAVPP"))
#' export_sq(sq_ami, "ape::AAbin")
#' export_sq(sq_ami, "Biostrings::AAStringSet", c("one", "two"))
#' export_sq(sq_ami, "seqinr::SeqFastaAA")
#' 
#' sq_nuc <- construct_sq(c("TGATGAAGCGCA", "TTGATGGGAA"))
#' export_sq(sq_nuc, "ape::DNAbin", name = c("one", "two"))
#' export_sq(sq_nuc, "Biostrings::DNAStringSet")
#' export_sq(sq_nuc, "seqinr::SeqFastadna")
#' @seealso \code{\link{sq}} \code{\link{import_sq}}
#' @export
export_sq <- function(sq, export_format, name) {
  .validate_sq(sq)
  if (!missing(name)) {
    .check_character(name, "'name'")
    .check_eq_lens(sq, name, "'sq'", "'name'")
  } else name <- NULL
  ami_formats <- c("seqinr::SeqFastaAA", "ape::AAbin", "Biostrings::AAStringSet")
  nuc_formats <- c("seqinr::SeqFastadna", "ape::DNAbin", "Biostrings::DNAStringSet")
  .check_export_format(export_format, ami_formats, nuc_formats)
  type <- .get_sq_type(sq)
  .check_type(type, "type of 'sq'")
  .check_type_matches_format(type, export_format, ami_formats, nuc_formats)
  
  if (export_format %in% c("seqinr::SeqFastaAA", "seqinr::SeqFastadna")) {
    .check_is_installed("seqinr")
    if (type == "ami") {
      ret <- lapply(.unpack_from_sq(sq, "char"), seqinr::as.SeqFastaAA)
    } else if (type %in% c("nuc", "dna", "rna")) {
      ret <- lapply(.unpack_from_sq(sq, "char"), seqinr::as.SeqFastadna)
    }
    
    if (!is.null(name)) {
      setNames(lapply(1:length(ret), function(i) {
        attr(ret[[i]], "name") <- name[i]
        ret[[i]]
      }), name)
    } else ret
  } else if (export_format %in% c("ape::AAbin", "ape::DNAbin")) {
    .check_is_installed("ape")
    if (type == "ami") {
      ape::as.AAbin(setNames(.unpack_from_sq(sq, "char"), name))
    } else if (type %in% c("nuc", "dna", "rna")) {
      ape::as.DNAbin(setNames(.unpack_from_sq(sq, "char"), name))
    }
  } else if (export_format %in% c("Biostrings::AAStringSet", "Biostrings::DNAStringSet")) {
    .check_is_installed("Biostrings")
    if (type == "ami") {
      Biostrings::AAStringSet(setNames(unlist(.unpack_from_sq(sq, "string")), name))
    } else if (type %in% c("nuc", "dna", "rna")) {
      Biostrings::DNAStringSet(setNames(unlist(.unpack_from_sq(sq, "string")), name))
    }
  }
}
