#' @exportMethod tbl_sum sqtbl
#' @export
tbl_sum.sqtbl <- function(x) {
  ret <- NextMethod()
  names(ret) <- "A sqtibble"
  ret
}

#' @importFrom pillar type_sum
#' @exportMethod type_sum sqcol
#' @export
type_sum.sqcol <- function(x) {
  "sq"
}

#' @exportMethod `[` sqcol
#' @export
`[.sqcol` <- function(i, ...) {
  ret <- NextMethod()
  class(ret) <- "sqcol"
  ret
}

#' @importFrom crayon col_nchar
#' @importFrom crayon green
#' @importFrom crayon silver
#' @importFrom pillar pillar_shaft
#' @importFrom pillar new_pillar_shaft
#' @importFrom pillar get_max_extent
#' @importFrom pillar get_extent
#' @exportMethod pillar_shaft sqcol
#' @export
pillar_shaft.sqcol <- function(x, ...) {
  x <- if (get_color_opt()) {
    as.character(lapply(x, function(i) paste(as.character(ifelse(is.na(i), silver("*"), green(as.character(i)))), collapse = "")))
  } else {
    as.character(lapply(x, function(i) paste(as.character(ifelse(is.na(i), "*", as.character(i))), collapse = "")))
  }
    
  longest_str <- get_max_extent(x)
  min_str_width <- ifelse(longest_str >= 6, 6, longest_str)
  
  opt <- get_print_length()
  
  new_pillar_shaft(x,
                   width = min(longest_str + col_nchar(longest_str) + 3, 
                               opt + col_nchar(longest_str) + 6),
                   min_width = col_nchar(longest_str) + min_str_width + 3,
                   class = "pillar_shaft_sqcol",
                   align = "left")
}

#' @importFrom crayon col_nchar
#' @importFrom crayon blue
#' @importFrom crayon silver
#' @importFrom crayon col_substring
#' @importFrom pillar new_ornament
#' @exportMethod format pillar_shaft_sqcol
#' @export
format.pillar_shaft_sqcol <- function(x, width, ...) {
  if (width < attr(x, "min_width")) {
    stop("need at least width ", attr(x, "min_width"), ", requested ", 
         width, ".", call. = FALSE)
  } 
  
  s_opt <- if (get_color_opt()) silver else identity
  b_opt <- if (get_color_opt()) blue else identity
  
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
  row <- paste0(row, whitespace_vec, b_opt(" ["), b_opt(seq_lengths), b_opt("]"))
  
  new_ornament(row, width = width, align = "left")
}

#' @exportMethod print sq
#' @export
print.sq <- function(x, ...) {
  sqtype <- intersect(class(x), c("aasq", "nucsq", "untsq"))
  if (length(sqtype) == 0 || length(sqtype) > 1) {
    sqtype <- "sq (improper subtype):\n"
  } else {
    sqtype <- paste0(c(aasq = "aa", nucsq = "nuc", untsq = "unt")[sqtype], " sequence:\n")
  }
  cat(paste0(sqtype, paste0(x, collapse = ""), "\n"))
}