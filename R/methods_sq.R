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
  ret <- do.call(rbind, lapply(x, function(row) row[1:max_len]))
  ret[ret == .get_na_char()] <- NA
  ret
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

#' @importFrom crayon blue
#' @importFrom crayon silver
#' @importFrom crayon green
#' @importFrom crayon col_nchar
#' @exportMethod print sq
#' @export
print.sq <- function(x,  
                     max_length = NULL,
                     max_sequences = getOption("tidysq_max_print_sequences"),
                     use_color = getOption("tidysq_colorful_sq_print"), 
                     letters_sep = NULL) {
  alph <- .get_alph(x)
  if (is.null(letters_sep)) {
    letters_sep <- if (all(nchar(alph) == 1)) "" else " "
  }
  p_width <- getOption("width")
  
  num_lines <- min(max_sequences, length(x))
  sq <- x[1:num_lines]
  sq_cut <- .cut_sq(sq, ceiling((p_width - 6) / (8 * (nchar(letters_sep) + 1))))
  sq_cut <- .debitify_sq(sq_cut, "char")
  
  na_char <- .get_na_char()
  if (use_color) sq_cut <- lapply(sq_cut, function(s) {s[s == na_char] <- silver(na_char); s})
  
  inds_width <- nchar(num_lines) + 2
  p_inds <- format(paste0("[", 1:num_lines, "]"),
              
                    width = inds_width, justify = "right")
  
  lens <- .get_lens(sq)
  lens_width <- max(nchar(lens)) + 2
  p_lens <- paste0("<", lens, ">")
  if (use_color) p_lens <- blue(p_lens)
  
  needs_dots <- rep(FALSE, num_lines)
  
  for (i in 1:num_lines) {
    s <- sq_cut[[i]]
    cum_lens <- cumsum(col_nchar(s)) + (0:(length(s) - 1)) * nchar(letters_sep)
    res_lens <- p_width - col_nchar(p_lens[i]) - nchar(p_inds[i]) - 2
    s <- s[cum_lens < res_lens]
    n <- length(s)
    if (n < lens[i]) {
      s <- s[cum_lens[1:n] < res_lens - 3]
      needs_dots[i] <- TRUE
    }
    sq_cut[[i]] <- s
  }
  
  p_body <- sapply(sq_cut, function(s) paste(s, collapse = letters_sep))
  if (use_color) p_body <- green(p_body)
  p_dots <- ifelse(needs_dots, "...", "")
  if (use_color) p_dots <- silver(p_dots)
  p_spaces <- sapply(p_width - col_nchar(p_body) - nchar(p_inds) - 
                       col_nchar(p_lens) - col_nchar(p_dots) - 2, 
                     function(l) paste0(rep(" ", l), collapse = ""))
  
  p_body <- paste0(p_inds, " ", p_body, p_dots, p_spaces, " ", p_lens, collapse = "\n")
  header <- .get_print_header(sq)
  footer <- .get_print_footer(x, num_lines)
  cat(header, p_body, footer, sep = "\n")
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

.cut_sq <- function(sq, num_oct) {
  alph_size <- tidysq:::.get_alph_size(tidysq:::.get_alph(sq))
  ret <- lapply(sq, function(s) if (length(s) < num_oct * alph_size) s 
         else s[1:(num_oct * alph_size)])
  .set_class_alph(ret, sq)
}

.get_print_header <- function(sq) {
  type <- .get_sq_type(sq)
  if (length(type) != 1) {
    "sq (improper subtype!):"
  } else {
    type_msg <- switch(type,
                       ami = "ami (amino acids)",
                       nuc = "nuc (nucleotides)",
                       unt = "unt (unspecified type)",
                       atp = "atp (atypical alphabet)")
    clean_msg <- if (.is_cleaned(sq)) ", cln (cleaned)" else ""
    paste0(type_msg, clean_msg, " sequences list:")
  }
}

.get_print_footer <- function(sq, num_lines) {
  if (length(sq) > num_lines) 
    paste0("printed ", num_lines, " out of ", length(sq), "")
  else ""
}