#methods other than print of class sq

#' @exportMethod `[` sq
#' @export
`[.sq` <- function(x, i, j, ...) {
  ret <- NextMethod()
  class(ret) <- class(x)
  attr(ret, "alphabet") <- .get_alph(x)
  ret
}

#' @exportMethod as.character sq
#' @export
as.character.sq <- function(x, ...) {
  alph <- .get_alph(x)
  sapply(.debitify_sq(x, alph), function(s) paste(ifelse(is.na(s), "*", s), collapse = ""))
}

#' @exportMethod is sq
#' @export
is.sq <- function(x) {
  tryCatch({validate_sq(x); TRUE}, error = function(e) FALSE)
}

#' @exportMethod is amisq
#' @export
is.amisq <- function(x) {
  tryCatch({validate_amisq(x); TRUE}, error = function(e) FALSE)
}

#' @exportMethod is nucsq
#' @export
is.nucsq <- function(x) {
  tryCatch({validate_nucsq(x); TRUE}, error = function(e) FALSE)
}

#' @exportMethod is untsq
#' @export
is.untsq <- function(x) {
  tryCatch({validate_untsq(x); TRUE}, error = function(e) FALSE)
}

#' @exportMethod is simsq
#' @export
is.simsq <- function(x) {
  tryCatch({validate_simsq(x); TRUE}, error = function(e) FALSE)
}

#' @exportMethod is atpsq
#' @export
is.atpsq <- function(x) {
  tryCatch({validate_atpsq(x); TRUE}, error = function(e) FALSE)
}
#' Compare sq object 
#' @description Compares input \code{\link{sq}} object with another given.
#'   
#' @details \code{`==``} checks if the input object is sequence, if yes, converts
#' it to chracters and checks whether given object can be compared with
#' \code{\link{sq}} object. If given sequence consists lowercase, the function
#' rewrites them into capital ones.
#' 
#' @param x1 \code{\link{sq}} object.
#' @param x2 an object to compare with \code{\link{sq}}.
#' 
#' @examples 
#' 
#' # Creating sq object to work on:
#' sq <- construct_sq(c("ACTGCTG", "CTTAGA", 
#'                      "CCCT", "CTGAATGT"), type = "nuc")
#' nuc_dna_sequence <- construct_sq(c("ACTGCTG", "CTTAGA", 
#'                                    "GGAA", "ATGAACGT"), type = "nuc")                              
#' # Get an overview of the sequences:
#' summary(sq)
#' summary(nuc_dna_sequence)
#' 
#' Comparing object nuc_dna_sequence to \code{\link{sq}}:
#' 
#' `==`(\code{\link{sq}},nuc_dna_sequence)
#'                                                                     
#' @exportMethod `==` sq
#' @export
`==.sq` <- function(e1, e2) {
  #TODO make it faster and lighter, maybe?
  if (is.sq(e2)) {
    e2 <- as.character(e2)
  } else if (!is.character(e2)) {
    stop ("you cannot compare 'sq' object to object that is not character vector or 'sq' object")
  }
  
  type <- .get_sq_type(e1)
  if (type %in% c("ami", "nuc")) {
    e2 <- toupper(e2)
  }
  
  as.character(e1) == e2
}

#' @exportMethod print sq
#' @export
print.sq <- function(x, ...) {
  sqclass <- .get_sq_subclass(x)
  cln_msg <- if (.is_cleaned(x)) " (cleaned)" else ""
  
  if (length(sqclass) != 1) {
    sqclass <- "sq (improper subtype!):\n"
  } else {
    sqclass <- paste0(c(amisq = "ami (amino acids)", 
                        nucsq = "nuc (nucleotides)", 
                        untsq = "unt (unspecified type)", 
                        simsq = "sim (simplified alphabet)",
                        atpsq = "atp (atypical alphabet)")[sqclass], cln_msg, " sequences vector:\n")
  }
  
  alph <- .get_alph(x)
  decoded <- .debitify_sq(x, alph)
  decoded <- sapply(decoded, function(s) ifelse(length(s) == 0, 
                                                "<NULL sq>", 
                                                paste(ifelse(is.na(s), "*", s), collapse = "")))
  max_width <- max(nchar(1:length(x)))
  inds <- paste0("[", 1:length(x), "] ")
  cat(sqclass, paste0(format(inds, width = max_width + 3, justify = "right"), 
                      decoded, 
                      collapse = "\n"), 
      "\n", sep = "")
}

#' @exportMethod print encsq
#' @export
print.encsq <- function(x, ...) {
  sqclass <- "enc (encoded values) sequences vector:\n"

  alph <- .get_alph(x)
  decoded <- .debitify_sq(x, alph)
  decoded <- sapply(decoded, function(s) ifelse(length(s) == 0, 
                                                "<NULL sq>", 
                                                paste(ifelse(is.na(s), "*", s), collapse = " ")))
  max_width <- max(nchar(1:length(x)))
  inds <- paste0("[", 1:length(x), "] ")
  cat(sqclass, paste0(format(inds, width = max_width + 3, justify = "right"), 
                      decoded, 
                      collapse = "\n"), 
      "\n", sep = "")
}
