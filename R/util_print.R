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
#' @exportMethod type_sum atpsq
#' @export
type_sum.atpsq <- function(x) {
  "atp"
}

#' @importFrom pillar type_sum
#' @exportMethod type_sum encsq
#' @export
type_sum.encsq <- function(x) {
  "enc"
}

#' @importFrom pillar type_sum
#' @exportMethod type_sum clnsq
#' @export
type_sum.clnsq <- function(x) {
  paste0("(c)", NextMethod())
}

#' Print sq object
#' 
#' @description Prints input \code{\link{sq}} object in a human-friendly form.  
#' 
#' @details \code{Print} method is used by default in each case of calling the 
#' \code{\link{sq}} object with default parameters. 
#' Only by explicit calling the \code{print} method parameters can be changed. 
#'  
#' \code{Print} checks if the input \code{\link{sq}} object is cleaned and includes 
#' this information alongside with type in the printed message. On the right side of 
#' the sequence, in angle brackets, the length of each sequence is printed (e.q. "<9>").
#' 
#' If the \code{max_sequences} parameter is supplied, the desired number of sequences 
#' is printed and this information is included in a message (e.q. "printed 1 out of 3"). 
#' Only \code{max_sequences} value smaller than the number of sequences in object 
#' affects the function. The default value indicating how many sequences should 
#' be printed is 10, but it can be changed in \code{\link[=tidysq-options]{package options}}. 
#' 
#' Default value of \code{use_color} parameter is \code{TRUE} - sequences are printed
#' in green and empty sequences, NA character and dots in gray. If this option is disabled, 
#' all sequences are in default color of console.
#' 
#' The \code{letters_sep} parameter indicates how the letters should be separated 
#' (they are not by default). Any character string can be supplied but 
#' \code{\link{NA_character_}}.
#' 
#' If sequences are too long, only leading characters are printed (as many as possible
#' in single line) and following dots indicating that sequence is truncated.
#' 
#' If sequences contain \code{\link{NA}} (‘Not Available’ / Missing Values) values, they 
#' are printed as "!" character, but it can be changed in 
#' \code{\link[=tidysq-options]{package options}}.
#' 
#' This is overloaded function from base package. It is selected when \code{\link{sq}} 
#' object is used as a parameter for print function. To see the generic function 
#' page, check \code{\link[base:print]{here}}.
#' 
#' @param x \code{\link{sq}} object
#' @param max_sequences \code{numeric} value indicating how many sequences 
#' should be printed
#' @param use_color \code{logical} value indicating if sequences should 
#' be colored
#' @param letters_sep \code{character} value indicating how the letters 
#' should be separated
#' @param ... 	further arguments passed to or from other methods. 
#' Unused.
#' 
#' @examples
#' 
#' # Creating sq objects using construct_sq:
#' sq_ami <- construct_sq(c("MIAANYTWIL","TIAALGNIIYRAIE", 
#'                          "NYERTGHLI", "MAYXXXIALN"), type = "ami")
#' sq_nuc <- construct_sq(c("ATGCAGGA", "GACCGAACGAN", 
#'                          "TGACGAGCTTA"), type = "nuc")
#' sq_unt <- construct_sq(c("ATGCAGGA!", "TGACGAGCTTA", "", "TIAALGNIIYRAIE"))
#' 
#' # Printing without explicit function calling with default parameters:
#' sq_ami
#' sq_nuc
#' sq_unt
#' 
#' # Printing with explicit function calling and specific parameters:
#' print(sq_ami)
#' print(sq_nuc, max_sequences = 1, use_color = FALSE)
#' print(sq_unt, letters_sep = ":")
#' 
#' # Printing of the cleaned object:
#' clean(sq_nuc)
#' print(clean(sq_nuc), letters_sep = "-", use_color = FALSE)
#' 
#' @seealso \code{\link{sq}} \code{\link{clean}} \code{\link{tidysq-options}}
#' 
#' @importFrom crayon blue
#' @importFrom crayon silver
#' @importFrom crayon green
#' @importFrom crayon col_nchar
#' @exportMethod print sq
#' @export
print.sq <- function(x,  
                     max_sequences = getOption("tidysq_max_print_sequences"),
                     use_color = getOption("tidysq_colorful_sq_print"), 
                     letters_sep = NULL, ...) {
  .check_integer(max_sequences, "'max_sequences'")
  .check_logical(use_color, "'use_color'")
  .check_character(letters_sep, "'letters_sep'", single_elem = TRUE, 
                   allow_zero_len = TRUE, allow_null = TRUE)
  alph <- .get_alph(x)
  
  #if parameter is NULL and all letters are length one, no space
  if (is.null(letters_sep)) {
    letters_sep <- if (all(nchar(alph) == 1)) "" else " "
  }
  p_width <- getOption("width")
  
  #select at most max_sequences to print
  num_lines <- min(max_sequences, length(x))
  if (num_lines == 0) {
    cat(.get_print_empty_sq(x))
    return()
  }  
  sq <- x[1:num_lines]
  
  #cut sq object so that we don't need to debitify long sequences
  # 6 is minimum length of p_lens and p_inds, 8 is byte length
  sq_cut <- .cut_sq(sq, ceiling((p_width - 6) / (8 * (nchar(letters_sep) + 1))))
  sq_cut <- .debitify_sq(sq_cut, "int")
  
  
  #color NA's
  na_char <- if (use_color) silver(.get_na_char()) else .get_na_char()
  na_val <- .get_na_val(alph)
  alph[na_val] <- na_char
  sq_cut <- lapply(sq_cut, function(s) {
    alph[s]
  })
  
  #max index number width
  inds_width <- nchar(num_lines) + 2
  
  #indices to print
  p_inds <- format(paste0("[", 1:num_lines, "]"), 
                   width = inds_width, justify = "right")
  
  #lengths of sequences
  lens <- .get_lens(sq)
  
  #max length number width
  lens_width <- max(nchar(lens)) + 2
  
  #lengths to print
  p_lens <- paste0("<", lens, ">")
  if (use_color) p_lens <- blue(p_lens)
  
  needs_dots <- rep(FALSE, num_lines)
  
  for (i in 1:num_lines) {
    if (lens[i] == 0) {
      sq_cut[[i]] <- "<NULL>"
    } else {
      s <- sq_cut[[i]]
      # we count how much characters can we print by counting cumulative extent
      cum_lens <- cumsum(col_nchar(s)) + (0:(length(s) - 1)) * nchar(letters_sep)
      
      #max length of this line is p_width minus lengths of lens and inds
      res_lens <- p_width - col_nchar(p_lens[i]) - nchar(p_inds[i]) - 2
      
      #we remove characters we cannot print
      s <- s[cum_lens < res_lens]
      n <- length(s)
      
      #if printed sequence is shorter than original, we need also space for dots
      if (n < lens[i]) {
        s <- s[cum_lens[1:n] < res_lens - 3]
        needs_dots[i] <- TRUE
      }
      sq_cut[[i]] <- s
    }
  }
  
  #paste sequene
  p_body <- sapply(sq_cut, function(s) paste(s, collapse = letters_sep))
  if (use_color) p_body <- sapply(1:num_lines, function(i) {
    if (lens[i] == 0) silver(p_body[i]) else green(p_body[i])
  })
  
  #dots
  p_dots <- ifelse(needs_dots, "...", "")
  if (use_color) p_dots <- silver(p_dots)
  
  #spaces between sequence and lens
  p_spaces <- sapply(p_width - col_nchar(p_body) - nchar(p_inds) - 
                       col_nchar(p_lens) - col_nchar(p_dots) - 2, 
                     function(l) paste0(rep(" ", l), collapse = ""))
  
  #paste and cat everything
  p_body <- paste0(p_inds, " ", p_body, p_dots, p_spaces, " ", p_lens, collapse = "\n")
  header <- .get_print_header(sq)
  footer <- .get_print_footer(x, num_lines)
  cat(header, p_body, footer, sep = "\n")
}

#' @importFrom crayon blue
#' @importFrom crayon silver
#' @importFrom crayon cyan
#' @importFrom crayon col_nchar
#' @exportMethod print encsq
#' @export
print.encsq <- function(x,
                        max_sequences = getOption("tidysq_max_print_sequences"),
                        use_color = getOption("tidysq_colorful_sq_print"), 
                        letters_sep = NULL,
                        digits = 2, ...) {
  .check_integer(max_sequences, "'max_sequences'")
  .check_logical(use_color, "'use_color'")
  .check_character(letters_sep, "'letters_sep'", single_elem = TRUE, 
                   allow_zero_len = TRUE, allow_null = TRUE)
  .check_integer(digits, "'digits'", allow_zero = TRUE)
  alph <- .get_alph(x)
  
  #if parameter is NULL default sep is space
  if (is.null(letters_sep)) {
    letters_sep <-  " "
  }
  p_width <- getOption("width")
  
  #select at most max_sequences to print
  num_lines <- min(max_sequences, length(x))
  sq <- x[1:num_lines]
  
  #format alphabet accordingly to passed params
  sq <- .set_alph(sq, format(alph, digits = digits, scientific = FALSE))
  
  #cut sq object so that we don't need to debitify long sequences
  # 6 is minimum length of p_lens and p_inds, 8 is byte length
   sq_cut <- .cut_sq(sq, ceiling((p_width - 6) / (8 * (nchar(letters_sep) + 1))))
  sq_cut <- .debitify_sq(sq_cut, "int")
  
  sq_cut <- lapply(sq_cut, function(s) {
    s <- alph[s]
    s[is.na(s)] <- "NA"
    s
  })
  
  
  #max index number width
  inds_width <- nchar(num_lines) + 2
  
  #indices to print
  p_inds <- format(paste0("[", 1:num_lines, "]"), 
                   width = inds_width, justify = "right")
  
  #lengths of sequences
  lens <- .get_lens(sq)
  
  #max length number width
  lens_width <- max(nchar(lens)) + 2
  
  #lengths to print
  p_lens <- paste0("<", lens, ">")
  if (use_color) p_lens <- blue(p_lens)
  
  needs_dots <- rep(FALSE, num_lines)
  
  for (i in 1:num_lines) {
    if (lens[i] == 0) {
      sq_cut[[i]] <- "<NULL>"
    } else {
      s <- sq_cut[[i]]
      # we count how much characters can we print by counting cumulative extent
      cum_lens <- cumsum(col_nchar(s)) + (0:(length(s) - 1)) * nchar(letters_sep)
      
      #max length of this line is p_width minus lengths of lens and inds
      res_lens <- p_width - col_nchar(p_lens[i]) - nchar(p_inds[i]) - 2
      
      #we remove characters we cannot print
      s <- s[cum_lens < res_lens]
      n <- length(s)
      
      #if printed sequence is shorter than original, we need also space for dots
      if (n < lens[i]) {
        s <- s[cum_lens[1:n] < res_lens - 3]
        needs_dots[i] <- TRUE
      }
      sq_cut[[i]] <- s
    }
  }
  
  #paste sequene
  p_body <- sapply(sq_cut, function(s) paste(s, collapse = letters_sep))
  if (use_color) p_body <- sapply(1:num_lines, function(i) {
    if (lens[i] == 0) silver(p_body[i]) else cyan(p_body[i])
  })
  
  #dots
  p_dots <- ifelse(needs_dots, "...", "")
  if (use_color) p_dots <- silver(p_dots)
  
  #spaces between sequence and lens
  p_spaces <- sapply(p_width - col_nchar(p_body) - nchar(p_inds) - 
                       col_nchar(p_lens) - col_nchar(p_dots) - 2, 
                     function(l) paste0(rep(" ", l), collapse = ""))
  
  #paste and cat everything
  p_body <- paste0(p_inds, " ", p_body, p_dots, p_spaces, " ", p_lens, collapse = "\n")
  header <- .get_print_header(sq)
  footer <- .get_print_footer(x, num_lines)
  cat(header, p_body, footer, sep = "\n")
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
  p_width <- getOption("width")
  letters_sep <- ""
  use_color <- .get_color_opt()
  alph <- .get_alph(x)
  
  #cut sq object so that we don't need to debitify long sequences
  # 6 is minimum length of p_lens and p_inds, 8 is byte length
  sq_cut <- .cut_sq(x, ceiling((p_width - 6) / (8 * (nchar(letters_sep) + 1))))
  sq_cut <- .debitify_sq(sq_cut, "int")
  
  #color NA's
  na_char <- if (use_color) silver(.get_na_char()) else .get_na_char()
  na_val <- .get_na_val(alph)
  alph[na_val] <- na_char
  sq_cut <- lapply(sq_cut, function(s) {
    alph[s]
  })
  
  lens <- .get_lens(x)
  max_len_width <- max(nchar(lens))
  
  max_str_width <- max(sapply(sq_cut, function(s) sum(col_nchar(s))))
  min_str_width <- if (max_str_width >= 6) 6 else max_str_width + 1
  
  opt <- .get_print_length()
  
  new_pillar_shaft(list(sq = sq_cut, lens = lens),
                   width = min(max_str_width + max_len_width + 4, 
                               opt + max_len_width + 7),
                   min_width = max_len_width + min_str_width + 3,
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
  
  lens <- x[["lens"]]
  x <- x[["sq"]]
  
  use_color <- .get_color_opt()
  
  #max length number width
  lens_width <- max(nchar(lens)) + 2
    
  #default letters_sep here is NULL
  letters_sep <- ""
  
  num_lines <- length(x)
  
  #lengths to print
  p_lens <- paste0("<", lens, ">")
  if (use_color) p_lens <- blue(p_lens)
  
  needs_dots <- rep(FALSE, num_lines)
  
  for (i in 1:num_lines) {
    if (lens[i] == 0) {
      x[[i]] <- "<NULL>"
    } else {
      s <- x[[i]]
      # we count how much characters can we print by counting cumulative extent
      cum_lens <- cumsum(col_nchar(s)) + (0:(length(s) - 1)) * nchar(letters_sep)
      
      #max length of this line is width minus lengths of lens 
      res_lens <- width - col_nchar(p_lens[i]) - 1
      
      #we remove characters we cannot print
      s <- s[cum_lens < res_lens]
      n <- length(s)
      
      #if printed sequence is shorter than original, we need also space for dots
      if (n < lens[i]) {
        s <- s[cum_lens[1:n] < res_lens - 3]
        needs_dots[i] <- TRUE
      }
      x[[i]] <- s
    }
  }
  
  #paste sequene
  p_body <- sapply(x, function(s) paste(s, collapse = letters_sep))
  if (use_color) p_body <- sapply(1:num_lines, function(i) {
    if (lens[i] == 0) silver(p_body[i]) else green(p_body[i])
  })
  
  #dots
  p_dots <- ifelse(needs_dots, "...", "")
  if (use_color) p_dots <- silver(p_dots)
  
  #spaces between sequence and lens
  p_spaces <- sapply(width - col_nchar(p_body) - 
                       col_nchar(p_lens) - col_nchar(p_dots) - 1, 
                     function(l) paste0(rep(" ", l), collapse = ""))
  
  #paste and cat everything
  p_body <- paste0(p_body, p_dots, p_spaces, " ", p_lens)
  
  new_ornament(p_body, width = width, align = "left")
}

#' @importFrom pillar pillar_shaft
#' @importFrom pillar new_pillar_shaft
#' @importFrom pillar get_max_extent
#' @exportMethod pillar_shaft encsq
#' @export
pillar_shaft.encsq <- function(x, ...) {
  p_width <- getOption("width")
  letters_sep <- " "
  use_color <- .get_color_opt()
  
  x <- .set_alph(x, format(.get_alph(x), digits = 1, scientific = FALSE))
  alph <- .get_alph(x)
  
  #cut sq object so that we don't need to debitify long sequences
  # 6 is minimum length of p_lens and p_inds, 8 is byte length
  sq_cut <- .cut_sq(x, ceiling((p_width - 6) / (8 * (nchar(letters_sep) + 1))))
  sq_cut <- .debitify_sq(sq_cut, "int")
  
  sq_cut <- lapply(sq_cut, function(s) {
    s <- alph[s]
    s[is.na(s)] <- "NA"
    s
  })
  
  lens <- .get_lens(x)
  max_len_width <- max(nchar(lens))
  
  max_str_width <- max(sapply(sq_cut, function(s) sum(col_nchar(s))))
  min_str_width <- if (max_str_width >= 6) 6 else max_str_width
  
  opt <- .get_print_length()
  
  new_pillar_shaft(list(sq = sq_cut, lens = lens),
                   width = min(max_str_width + max_len_width + 3, 
                               opt + max_len_width + 6),
                   min_width = max_len_width + min_str_width + 3,
                   class = "pillar_shaft_encsq",
                   align = "left")
}

#' @importFrom crayon cyan
#' @importFrom crayon silver
#' @importFrom pillar new_ornament
#' @exportMethod format pillar_shaft_encsq
#' @export
format.pillar_shaft_encsq <- function(x, width, ...) {
  if (width < attr(x, "min_width")) {
    stop("need at least width ", attr(x, "min_width"), ", requested ", 
         width, ".", call. = FALSE)
  } 
  
  lens <- x[["lens"]]
  x <- x[["sq"]]
  
  use_color <- .get_color_opt()
  
  #max length number width
  lens_width <- max(nchar(lens)) + 2
    
  #default letters_sep here is NULL
  letters_sep <- " "
  
  num_lines <- length(x)
  
  #lengths to print
  p_lens <- paste0("<", lens, ">")
  if (use_color) p_lens <- blue(p_lens)
  
  needs_dots <- rep(FALSE, num_lines)
  
  for (i in 1:num_lines) {
    if (lens[i] == 0) {
      x[[i]] <- "<NULL>"
    } else {
      s <- x[[i]]
      # we count how much characters can we print by counting cumulative extent
      cum_lens <- cumsum(col_nchar(s)) + (0:(length(s) - 1)) * nchar(letters_sep)
      
      #max length of this line is width minus lengths of lens 
      res_lens <- width - col_nchar(p_lens[i]) - 1
      
      #we remove characters we cannot print
      s <- s[cum_lens < res_lens]
      n <- length(s)
      
      #if printed sequence is shorter than original, we need also space for dots
      if (n < lens[i]) {
        s <- s[cum_lens[1:n] < res_lens - 3]
        needs_dots[i] <- TRUE
      }
      x[[i]] <- s
    }
  }
  
  #paste sequene
  p_body <- sapply(x, function(s) paste(s, collapse = letters_sep))
  if (use_color) p_body <- sapply(1:num_lines, function(i) {
    if (lens[i] == 0) silver(p_body[i]) else cyan(p_body[i])
  })
  
  #dots
  p_dots <- ifelse(needs_dots, "...", "")
  if (use_color) p_dots <- silver(p_dots)
  
  #spaces between sequence and lens
  p_spaces <- sapply(width - col_nchar(p_body) - 
                       col_nchar(p_lens) - col_nchar(p_dots) - 1, 
                     function(l) paste0(rep(" ", l), collapse = ""))
  
  #paste and cat everything
  p_body <- paste0(p_body, p_dots, p_spaces, " ", p_lens)
  
  new_ornament(p_body, width = width, align = "left")
}

.cut_sq <- function(sq, num_oct) {
  alph_size <- .get_alph_size(.get_alph(sq))
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
                       atp = "atp (atypical alphabet)",
                       enc = "enc (encoded values)")
    clean_msg <- if (.is_cleaned(sq)) ", cln (cleaned)" else ""
    paste0(type_msg, clean_msg, " sequences list:")
  }
}

.get_print_footer <- function(sq, num_lines) {
  if (length(sq) > num_lines) 
    paste0("printed ", num_lines, " out of ", length(sq), "")
  else ""
}

.get_print_empty_sq <- function(sq) {
  type <- .get_sq_type(sq)
  if (length(type) != 1) {
    "sq (improper subtype!):"
  } else {
    type_msg <- switch(type,
                       ami = "ami (amino acids)",
                       nuc = "nuc (nucleotides)",
                       unt = "unt (unspecified type)",
                       atp = "atp (atypical alphabet)",
                       enc = "enc (encoded values)")
    clean_msg <- if (.is_cleaned(sq)) ", cln (cleaned)" else ""
    paste0(type_msg, clean_msg, " sequences list of length 0")
  }
}