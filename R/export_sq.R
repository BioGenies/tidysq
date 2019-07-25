#' @export
export_sq <- function(sq, export_format, name) {
  validate_sq(sq)
  if (!missing(name) && 
      (!is.character(name) ||
       any(is.null(name) | is.na(name)) ||
       (length(name) != length(sq)))) {
    stop("'name' has to be a non-NULL character vector without NA's, of lenght equal to length of sq")
  }
  if (missing(name)) name = NULL
  
  ami_formats <- c("seqinr::SeqFastaAA", "ape::AAbin", "Biostrings::AAStringSet")
  nuc_formats <- c("seqinr::SeqFastadna", "ape::DNAbin", "Biostrings::DNAStringSet")
  if (missing(export_format) ||
      !(export_format %in% c(ami_formats, nuc_formats))) {
    stop("you need to specify proper 'export_format'; check manual for possible formats")
  }
  type <- .get_sq_type(sq)
  if (!(type %in% c("ami", "nuc"))) {
    stop("'sq' has to have 'ami' or 'nuc' type")
  }
  if (type == "ami" && !(export_format %in% ami_formats) ||
      type == "nuc" && !(export_format %in% nuc_formats)) {
    stop("'sq' object type doesn't match 'export_format'")
  }
  
  if (export_format %in% c("seqinr::SeqFastaAA", "seqinr::SeqFastadna")) {
    if (!("seqinr" %in% rownames(installed.packages()))) {
      stop("you need to install 'seqinr' package to export object to its formats")
    }
    if (type == "ami") {
      ret <- lapply(.debitify_sq(sq, .get_alph(sq)), seqinr::as.SeqFastaAA)
    } else if (type == "nuc") {
      ret <- lapply(.debitify_sq(sq, .get_alph(sq)), seqinr::as.SeqFastadna)
    }
    
    if (!is.null(name)) {
      setNames(lapply(1:length(ret), function(i) {
        attr(ret[[i]], "name") <- name[i]
        ret[[i]]
      }), name)
    } else ret
  } else if (export_format %in% c("ape::AAbin", "ape::DNAbin")) {
    if (!("ape" %in% rownames(installed.packages()))) {
      stop("you need to install 'ape' package to export object to its formats")
    } 
    if (type == "ami") {
      ape::as.AAbin(setNames(.debitify_sq(sq, .get_alph(sq)), name))
    } else if (type == "nuc") {
      ape::as.DNAbin(setNames(.debitify_sq(sq, .get_alph(sq)), name))
    }
  } else if (export_format %in% c("Biostrings::AAStringSet", "Biostrings::DNAStringSet")) {
    if (!("Biostrings" %in% rownames(installed.packages()))) {
      stop("you need to install 'Biostrings' package to export object to its formats")
    } 
    if (type == "ami") {
      Biostrings::AAStringSet(setNames(sq, name))
    } else if (type == "nuc") {
      Biostrings::DNAStringSet(setNames(sq, name))
    }
  }
}