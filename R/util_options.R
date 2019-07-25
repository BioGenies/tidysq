#internal functions for accessing options of package
.handle_opt_txt <- function(option_name, txt) {
  opt <- getOption(option_name)
  
  opt <- ifelse(is.null(opt) || 
            !is.character(opt) || 
            length(opt) != 1 || 
            !(opt %in% c("error", "warning", "message", "none")), "warning", opt)
  
  switch(opt,
         error = stop(txt, call. = FALSE),
         warning = warning(txt, call. = FALSE),
         message = message(txt),
         none = invisible())
}

.get_print_length <- function() {
  opt <- getOption("tidysq_max_sq_print_width")
  
  ifelse (is.null(opt) || 
            !is.numeric(opt) || 
            is.na(opt) || 
            is.nan(opt) || 
            is.infinite(opt) || 
            !(opt > 0), 15, opt)
}

.get_color_opt <- function() {
  opt <- getOption("tidysq_colorful_sq_print")
  
  ifelse (is.null(opt) || 
            !is.logical(opt) || 
            is.na(opt), FALSE, opt)
}

.get_na_char <- function() {
  opt <- getOption("tidysq_na_print_char")

  ifelse (is.null(opt) ||
            !is.character(opt) ||
            is.na(opt) ||
            (length(opt) != 1), "!", opt)  
}