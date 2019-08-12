#' @export
export_sq <- function(sq, export_format, name) {
  validate_sq(sq)
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
      ret <- lapply(.debitify_sq(sq, "char"), seqinr::as.SeqFastaAA)
    } else if (type == "nuc") {
      ret <- lapply(.debitify_sq(sq, "char"), seqinr::as.SeqFastadna)
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
      ape::as.AAbin(setNames(.debitify_sq(sq, "char"), name))
    } else if (type == "nuc") {
      ape::as.DNAbin(setNames(.debitify_sq(sq, "char"), name))
    }
  } else if (export_format %in% c("Biostrings::AAStringSet", "Biostrings::DNAStringSet")) {
    .check_is_installed("Biostrings")
    if (type == "ami") {
      Biostrings::AAStringSet(setNames(unlist(.debitify_sq(sq, "string")), name))
    } else if (type == "nuc") {
      Biostrings::DNAStringSet(setNames(unlist(.debitify_sq(sq, "string")), name))
    }
  }
}
