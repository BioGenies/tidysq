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
  .debitify_sq(x, "string")
}

#' @exportMethod as.matrix sq
#' @export
as.matrix.sq <- function(x, ...) {
  x <- .debitify_sq(x, "char")
  max_len <- max(lengths(x))
  do.call(rbind, lapply(x, function(row) row[1:max_len]))
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

#' @exportMethod print sq
#' @export
print.sq <- function(x, ...) {
  sqclass <- .get_sq_subclass(x)
  cln_msg <- if (.is_cleaned(x)) " (cleaned)" else ""
  na_char <- .get_na_char()
  
  if (length(sqclass) != 1) {
    sqclass <- "sq (improper subtype!):\n"
  } else {
    sqclass <- paste0(c(amisq = "ami (amino acids)", 
                        nucsq = "nuc (nucleotides)", 
                        untsq = "unt (unspecified type)", 
                        atpsq = "atp (atypical alphabet)")[sqclass], cln_msg, " sequences vector:\n")
  }
  
  alph <- .get_alph(x)
  decoded <- .debitify_sq(x, "string")
  decoded <- sapply(decoded, function(s) ifelse(s == "" , "<NULL sq>", s))
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
  na_char <- .get_na_char()

  alph <- .get_alph(x)
  decoded <- .apply_sq(x, "int", "none", function(s) alph[s])
  decoded <- sapply(decoded, function(s) ifelse(length(s) == 0, 
                                                "<NULL sq>", 
                                                paste(ifelse(is.na(s), na_char, s), collapse = " ")))
  max_width <- max(nchar(1:length(x)))
  inds <- paste0("[", 1:length(x), "] ")
  cat(sqclass, paste0(format(inds, width = max_width + 3, justify = "right"), 
                      decoded, 
                      collapse = "\n"), 
      "\n", sep = "")
}
