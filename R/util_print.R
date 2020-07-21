#' @importFrom pillar type_sum
#' @export
type_sum.amisq <- function(x) {
  "ami"
}

#' @importFrom pillar type_sum
#' @export
type_sum.dnasq <- function(x) {
  "dna"
}

#' @importFrom pillar type_sum
#' @export
type_sum.rnasq <- function(x) {
  "rna"
}

#' @importFrom pillar type_sum
#' @export
type_sum.untsq <- function(x) {
  "unt"
}

#' @importFrom pillar type_sum
#' @export
type_sum.atpsq <- function(x) {
  "atp"
}

#' @importFrom pillar type_sum
#' @export
type_sum.encsq <- function(x) {
  "enc"
}

#' @importFrom pillar type_sum
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
#' in green, while empty sequences, NA character and dots in gray. If this option is disabled, 
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
#' sq_dna <- construct_sq(c("ATGCAGGA", "GACCGAACGAN", 
#'                          "TGACGAGCTTA"), type = "dna")
#' sq_unt <- construct_sq(c("ATGCAGGA!", "TGACGAGCTTA", "", "TIAALGNIIYRAIE"))
#' 
#' # Printing without explicit function calling with default parameters:
#' sq_ami
#' sq_dna
#' sq_unt
#' 
#' # Printing with explicit function calling and specific parameters:
#' print(sq_ami)
#' print(sq_dna, max_sequences = 1, use_color = FALSE)
#' print(sq_unt, letters_sep = ":")
#' 
#' # Printing of the cleaned object:
#' clean(sq_dna)
#' print(clean(sq_dna), letters_sep = "-", use_color = FALSE)
#' 
#' @seealso \code{\link{sq}} \code{\link{clean}} \code{\link{tidysq-options}}
#' 
#' @importFrom crayon blue
#' @importFrom crayon silver
#' @importFrom crayon green
#' @importFrom crayon col_nchar
#' @export
print.sq <- function(x,  
                     max_sequences = getOption("tidysq_p_max_sequences"),
                     use_color = getOption("tidysq_p_use_color"),
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
  # TODO: maybe replace it with "char" or "string"? that is, considering the need of coloring NA's
  sq_cut <- .unpack_from_sq(sq_cut, "int")
  
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
  lens <- lengths(sq)
  
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
  
  #paste sequence
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
  header <- .get_print_nonempty_header(sq)
  footer <- .get_print_footer(x, num_lines)
  cat(header, p_body, if (length(x) > num_lines) footer, sep = "\n")
}

#' @importFrom crayon blue
#' @importFrom crayon silver
#' @importFrom crayon cyan
#' @importFrom crayon col_nchar
#' @export
print.encsq <- function(x,
                        max_sequences = getOption("tidysq_p_max_sequences"),
                        use_color = getOption("tidysq_p_use_color"),
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
  sq_cut <- .unpack_from_sq(sq_cut, "int")
  
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
  lens <- lengths(sq)
  
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
  
  #paste sequence
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
  header <- .get_print_nonempty_header(sq)
  footer <- .get_print_footer(x, num_lines)
  cat(header, p_body, if (length(x) > num_lines) footer, sep = "\n")
}

#' @importFrom crayon green
#' @importFrom crayon silver
#' @importFrom pillar pillar_shaft
#' @export
pillar_shaft.sq <- function(x, ...) {
  # color NA's
  alph <- .get_alph(x)
  na_char <- if (.get_color_opt()) silver(.get_na_char()) else .get_na_char()
  na_val <- .get_na_val(alph)
  alph[na_val] <- na_char
  
  .pillar_shaft_sq(x, alph, "", green)
}

#' @importFrom crayon cyan
#' @importFrom pillar pillar_shaft
#' @export
pillar_shaft.encsq <- function(x, ...) {
  x <- .set_alph(x, format(.get_alph(x), digits = 1, scientific = FALSE))
  alph <- .get_alph(x)
  
  .pillar_shaft_sq(x, alph, " ", cyan)
}

#' @importFrom crayon col_nchar
#' @importFrom crayon blue
#' @importFrom crayon silver
#' @importFrom pillar new_ornament
#' @export
format.pillar_shaft_sq <- function(x, width, ...) {
  if (width < attr(x, "min_width"))
    stop("need at least width ", attr(x, "min_width"), ", requested ", width, ".", call. = FALSE)
  
  lens <- attr(x, "lens")
  letters_sep <- attr(x, "letters_sep")
  body_color <- attr(x, "body_color")
  align <- attr(x, "align")
  
  # max width of length number
  lens_width <- max(nchar(lens)) + 2
  
  x <- mapply(function(sq, len) {
    if (len == 0) {
      structure("<NULL>", dots = "")
    } else {
      # we count how much characters can we print by counting cumulative extent
      cum_lens <- cumsum(col_nchar(sq)) + (0:(length(sq) - 1)) * nchar(letters_sep)
      # max length of this line is its width minus the lens_width
      res_lens <- width - lens_width - 1
      
      # if total length is greater than reserved space, we have to cut it down
      if (tail(cum_lens, n = 1) > res_lens) {
        # if printed sequence is shorter than original, we also need space for dots
        # find last index that allows sq to fit into reserved space
        n <- Position(identity, cum_lens <= res_lens - 3, right = TRUE)
        sq <- sq[1:n]
        attr(sq, "dots") <- "..."
      } else {
        attr(sq, "dots") <- ""
      }
      sq
    }
  }, x, lens, SIMPLIFY = FALSE)
  
  # strings formatted (without colors) to be printed
  # sequences
  p_body <- sapply(x, paste, collapse = letters_sep)
  # dots
  p_dots <- sapply(x, attr, "dots")
  # lengths
  p_lens <- paste0("<", lens, ">")
  # spaces between sequences and lens
  p_spaces <- sapply(width - col_nchar(p_body) - col_nchar(p_lens) - col_nchar(p_dots),
                     function(l) paste0(rep(" ", l), collapse = ""))
  
  # add colors to text if necessary
  if (.get_color_opt()) {
    p_lens <- blue(p_lens)
    p_body <- mapply(function(content, len)
      if (len == 0) silver(content) else body_color(content),
      p_body, lens)
    p_dots <- silver(p_dots)
  }
  
  # paste and cat everything
  new_ornament(paste0(p_body, p_dots, p_spaces, p_lens),
               width = width, align = align)
}

.cut_sq <- function(sq, num_oct) {
  # note: num_oct should be greater than 0, not that it'd make sense to use 0 or less
  alph_size <- .get_alph_size(.get_alph(sq))
  ret <- lapply(sq, function(s) {
    if (length(s) <= num_oct * alph_size) s 
    else structure(s[1:(num_oct * alph_size)], original_length = attr(s, "original_length"))
  })
  .set_class_alph(ret, sq)
}

.get_print_header <- function(sq) {
  type <- .get_sq_type(sq)
  if (length(type) != 1) {
    "sq (improper subtype!)"
  } else {
    type_msg <- switch(type,
                       ami = "ami (amino acids)",
                       dna = "dna (DNA)",
                       rna = "rna (RNA)",
                       unt = "unt (unspecified type)",
                       atp = "atp (atypical alphabet)",
                       enc = "enc (encoded values)")
    clean_msg <- if (.is_cleaned(sq)) ", cln (cleaned)" else ""
    paste0(type_msg, clean_msg, " sequences list")
  }
}

.get_print_nonempty_header <- function(sq) {
  paste0(.get_print_header(sq), ":")
}

.get_print_empty_sq <- function(sq) {
  paste0(.get_print_header(sq), " of length 0")
}

.get_print_footer <- function(sq, num_lines) {
  paste0("printed ", num_lines, " out of ", length(sq), "")
}

#' @importFrom crayon col_nchar
#' @importFrom pillar new_pillar_shaft
.pillar_shaft_sq <- function(x, alph, letters_sep, body_color) {
  p_width <- getOption("width")
  
  # cut sq object so that we don't need to debitify long sequences
  # 6 is minimum length of p_lens and p_inds, 8 is byte length
  sq_cut <- .cut_sq(x, ceiling((p_width - 6) / (8 * (nchar(letters_sep) + 1))))
  sq_cut <- .unpack_from_sq(sq_cut, "int")
  sq_cut <- lapply(sq_cut, function(s) alph[s])
  
  # maximum length of length numbers
  lens <- lengths(x)
  max_len_width <- max(nchar(lens))
  
  max_str_width <- max(sapply(sq_cut, function(s) sum(col_nchar(s))))
  min_str_width <- if (max_str_width >= 6) 6 else max_str_width
  
  new_pillar_shaft(sq_cut,
                   width = max_len_width + 3 + min(max_str_width, .get_print_length() + 3),
                   min_width = max_len_width + 3 + min_str_width,
                   class = "pillar_shaft_sq",
                   align = "left",
                   lens = lens,
                   letters_sep = letters_sep,
                   body_color = body_color)
}
