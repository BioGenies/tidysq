#' @exportMethod tbl_sum sqtbl
#' @export
tbl_sum.sqtbl <- function(x) {
  ret <- NextMethod()
  names(ret) <- "A sqtibble"
  ret
}

#' @importFrom pillar type_sum
#' @exportMethod type_sum amisq
#' @export
type_sum.amisq <- function(x) {
  "ami"
}

#' @importFrom pillar type_sum
#' @exportMethod type_sum nucsq
#' @export
type_sum.nucsq <- function(x) {
  "nuc"
}

#' @importFrom pillar type_sum
#' @exportMethod type_sum untsq
#' @export
type_sum.untsq <- function(x) {
  "unt"
}

#' @importFrom pillar type_sum
#' @exportMethod type_sum simsq
#' @export
type_sum.simsq <- function(x) {
  "sim"
}

#' @importFrom pillar type_sum
#' @exportMethod type_sum atpsq
#' @export
type_sum.atpsq <- function(x) {
  "atp"
}

#' @importFrom pillar type_sum
#' @exportMethod type_sum clnsq
#' @export
type_sum.clnsq <- function(x) {
  paste0("(c)", NextMethod())
}

#' @importFrom crayon col_nchar
#' @importFrom crayon green
#' @importFrom crayon silver
#' @importFrom pillar pillar_shaft
#' @importFrom pillar new_pillar_shaft
#' @importFrom pillar get_max_extent
#' @importFrom pillar get_extent
#' @exportMethod pillar_shaft sq
#' @export
pillar_shaft.sq <- function(x, ...) {
  alph <- .get_alph(x)
  x <- lapply(x, function(s) alph[s])
  x <- if (.get_color_opt()) {
    lapply(x, function(s) paste(
      as.character(ifelse(is.na(s), 
                          silver("*"), 
                          green(s))), 
      collapse = ""))
  } else {
    lapply(x, function(s) paste(
      as.character(ifelse(is.na(s), 
                          "*", 
                          s)), 
      collapse = ""))
  }
    
  longest_str <- get_max_extent(x)
  min_str_width <- if (longest_str >= 6) 6 else longest_str
  
  opt <- .get_print_length()
  
  new_pillar_shaft(x,
                   width = min(longest_str + col_nchar(longest_str) + 3, 
                               opt + col_nchar(longest_str) + 6),
                   min_width = col_nchar(longest_str) + min_str_width + 3,
                   class = "pillar_shaft_sq",
                   align = "left")
}

#' @importFrom crayon col_nchar
#' @importFrom crayon blue
#' @importFrom crayon silver
#' @importFrom crayon col_substring
#' @importFrom pillar new_ornament
#' @exportMethod format pillar_shaft_sq
#' @export
format.pillar_shaft_sq <- function(x, width, ...) {
  if (width < attr(x, "min_width")) {
    stop("need at least width ", attr(x, "min_width"), ", requested ", 
         width, ".", call. = FALSE)
  } 
  
  s_opt <- if (.get_color_opt()) silver else identity
  b_opt <- if (.get_color_opt()) blue else identity
  
  seq_lengths <- col_nchar(x)
  nums_lengths <- nchar(seq_lengths)
  available_seq_space <- width - nums_lengths - 3
  inds_need_dots <- available_seq_space < seq_lengths
  available_seq_space[inds_need_dots] <- available_seq_space[inds_need_dots] - 3
  row <- col_substring(x, 1, available_seq_space)
  row[inds_need_dots] <- paste0(row[inds_need_dots], s_opt("..."))
  whitespace_lenghts <- available_seq_space - seq_lengths
  whitespace_lenghts[whitespace_lenghts < 0] <- 0
  whitespace_vec <- sapply(whitespace_lenghts, 
                           function(n) paste(rep(" ", n), collapse = ""))
  row <- paste0(row, whitespace_vec, b_opt(" <"), b_opt(seq_lengths), b_opt(">"))
  
  new_ornament(row, width = width, align = "left")
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
  
  dict <- .get_alph(x)
  names(dict) <- 1:length(dict)
  decoded <- sapply(x, function(s) paste(dict[s], collapse = ""))
  cat(paste0(sqclass, paste0("[", 1:length(x), "]  ", decoded, collapse = "\n")))
}