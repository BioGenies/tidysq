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

#' @importFrom pillar pillar_shaft
#' @importFrom pillar new_pillar_shaft
#' @importFrom pillar get_max_extent
#' @importFrom pillar get_extent
#' @exportMethod pillar_shaft sqcol
#' @export
pillar_shaft.sqcol <- function(x, ...) {
  x <- as.character(sapply(x, function(i) paste(as.character(ifelse(is.na(i), "*", as.character(i))), collapse = "")))
  # # x <- format(x, formatter = sqcol_formatter)
  # inds_too_long <- (nchar(x) > 20)
  # x[inds_too_long] <- paste0(strtrim(x[inds_too_long], 17), "...")
  longest_str <- get_max_extent(x)
  min_str_width <- ifelse(longest_str >= 6, 6, longest_str)
  
  new_pillar_shaft(x,
                   width = longest_str + nchar(longest_str) + 3,
                   min_width = nchar(longest_str) + min_str_width + 3,
                   class = "pillar_shaft_sqcol",
                   align = "left")
}

#' @importFrom pillar new_ornament
#' @exportMethod format pillar_shaft_sqcol
#' @export
format.pillar_shaft_sqcol <- function(x, width, ...) {
  if (width < attr(x, "min_width")) {
    stop("need at least width ", attr(x, "min_width"), ", requested ", 
         width, ".", call. = FALSE)
  }
  
  seq_lengths <- nchar(x)
  nums_lengths <- nchar(seq_lengths)
  available_seq_space <- width - nums_lengths - 3
  inds_need_dots <- available_seq_space < seq_lengths
  available_seq_space[inds_need_dots] <- available_seq_space[inds_need_dots] - 3
  row <- strtrim(x, available_seq_space)
  row[inds_need_dots] <- paste0(row[inds_need_dots], "...")
  whitespace_lenghts <- available_seq_space - seq_lengths
  whitespace_lenghts[whitespace_lenghts < 0] <- 0
  whitespace_vec <- sapply(whitespace_lenghts, 
                           function(n) paste(rep(" ", n), collapse = ""))
  row <- paste0(row, whitespace_vec," [", seq_lengths, "]")
  
  new_ornament(row, width = width, align = "left")
}

