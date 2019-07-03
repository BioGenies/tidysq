#' @exportClass nucsq
#' @export
construct_nucsq <- function(sq) {
  stopifnot(is.character(sq))
  sq <- toupper(sq)
  if(length(sq) == 1) {
    sq <- strsplit(sq, "")[[1]]
  } else {
    stopifnot(all(nchar(sq) == 1))
  }
  is_nuc_sq <- all(sq %in% nucleotides_df[,"one"])
  stopifnot(is_nuc_sq)
  
  object <- factor(sq, levels = nucleotides_df[,"one"])
  class(object) <- c("nuc_sq", "sq", class(object))
  object
}

#' @exportClass aasq
#' @export
construct_aasq <- function(sq) {
  # TO DO: what if user gives list of aminoacids three letter names
  stopifnot(is.character(sq))
  sq <- toupper(sq)
  if(length(sq) == 1) {
    sq <- strsplit(sq, "")[[1]]
  } else {
    stopifnot(all(nchar(sq) == 1))
  }
  is_aa_sq <- all(sq %in% aminoacids_df[,"one"])
  stopifnot(is_aa_sq)
  
  object <- factor(sq, levels = aminoacids_df[,"one"])
  class(object) <- c("aa_sq", "sq", class(object))
  object
}

#' @exportClass nucsq
#' @export
construct_nucsq <- function(sq) {
  stopifnot(is.character(sq))
  sq <- toupper(sq)
  if(length(sq) == 1) {
    sq <- strsplit(sq, "")[[1]]
  } else {
    stopifnot(all(nchar(sq) == 1))
  }
  is_nuc_sq <- all(sq %in% nucleotides_df[,"one"])
  stopifnot(is_nuc_sq)
  
  object <- factor(sq, levels = nucleotides_df[,"one"])
  class(object) <- c("nuc_sq", "sq", class(object))
  object
}