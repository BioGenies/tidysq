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

#' Print sq object
#' 
#' @description Prints input \code{\link{sq}} object in a human-friendly form.  
#' 
#' @details \code{Print} checks if the input \code{\link{sq}} object is cleaned and includes this information alongside with type in the printed message. 
#' All \code{\link{NA}} values are replaced with '*' and empty sequences are distinguished.
#' 
#' \code{Print} method is used by default in each case of calling the \code{\link{sq}} object.
#' 
#' This is overloaded function from base package. It is selected when \code{\link{sq}} object is used as a parameter for print function. To see the generic function page, 
#' check \link[here]{https://www.rdocumentation.org/packages/base/versions/3.6.1/topics/print}.
#' 
#' @param x \code{\link{sq}} object.
#' @param ... further arguments passed to or from other methods.
#' 
#' @examples
#' 
#' # Creating sq object to work on:
#' sq <- construct_sq(c("fafasfasfFSA", "ygagayagfa", "adsDaf"), type = "ami")
#' 
#' # Printing without implicit function calling:
#' sq
#' 
#' # Printing with implicit function calling:
#' print(sq)
#' 
#' # Implicit printing of the uncleaned object:
#' print(construct_sq("ACTAGAGTGATAGTAGGAGTAGA", type = "nuc"))
#'
#' # Implicit printing of the cleaned object:
#' print(clean(construct_sq("ACTAGAGTGATAGTAGGAGTAGA", type = "nuc")))
#' 
#' # Implicit printing of the object without defined type:
#' print(construct_sq(c("afsfd", "q243faadfa", "afsw34gesfv", "adfq2", "fasfas", "g'qp9u2r3'b;")))
#' 
#' # Implicit printing of the object with empty sequence:
#' print(construct_sq(c("afsfd", "", "adfq2", "fasfas", "")))
#' 
#' # Implicit printing of the object with NA element:
#' print(construct_sq(c("afsfd", NA, "adfq2", NA, "")))
#' 
#' # Implicit printing of the simplified object:
#' enc <- c(A = "a", B = "a", C = "a", D = "a", E = "a", F = "b", G = "b", 
#'          H = "b", I = "c", J = "c", K = "c", L = "c", M = "c", N = "c", 
#'          O = "c", P = "d", Q = "d", R = "d", S = "d", T = "d", U = "d", 
#' V = "d", W = "d", X = "d", Y = "d", Z = "d", `-` = "d")
#' print(simplify(sq, enc))
#'  
#' @seealso \link{sq} \link{clean} \link{simplify}
#' 
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
