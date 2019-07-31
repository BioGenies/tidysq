#' @export
export_sq <- function(sq, export_format, name) {
  validate_sq(sq)
  if (!missing(name) && 
      (!is.character(name) ||
       any(is.null(name) | is.na(name)) ||
       (length(name) != length(sq)))) {
    stop("'name' has to be a non-NULL character vector without NA's, of length equal to length of sq")
  }
  ami_formats <- c("seqinr::SeqFastaAA", "ape::AAbin", "Bioconductor::AAStringSet")
  nuc_formats <- c("seqinr::SeqFastadna", "ape::DNAbin", "Bioconductor::DNAStringSet")
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
    if (!("seqinr" %in% installed.packages()[["Package"]])) {
      stop("you need to install 'ape' package to export object to its formats")
    }
    if (type == "ami") {
      seqinr::as.SeqFastaAA(strsplit(as.character(sq)), name = name)
    }
    if (type == "nuc") {
      seqinr::as.SeqFastadna(strsplit(as.character(sq)), name = name)
    }
  }
  
  if (export_format %in% c("ape::AAbin", "ape::DNAbin")) {
    if (!("ape" %in% installed.packages()[["Package"]])) {
      stop("you need to install 'ape' package to export object to its formats")
    } 
    if (type == "ami") {
      ape::as.AAbin(setNames(as.list(as.character(sq)), name))
    }
    if (type == "nuc") {
      ape::as.DNAbin(setNames(as.list(as.character(sq)), name))
    }
  }
  
  if (export_format %in% c("ape::AAbin", "ape::DNAbin")) {
    if (!("Biostrings" %in% installed.packages()[["Package"]])) {
      stop("you need to install 'ape' package to export object to its formats")
    } 
    if (type == "ami") {
      Biostrings::AAStringSet(setNames(sq, name))
    }
    if (type == "nuc") {
      Biostrings::DNAStringSet(setNames(sq, name))
    }
  }
  
}
