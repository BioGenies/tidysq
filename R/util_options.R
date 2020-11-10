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

.is_fast_mode <- function() {
  opt <- getOption("tidysq_g_fast_mode")
  
  ifelse (is.null(opt) || 
            !is.logical(opt) || 
            is.na(opt), FALSE, opt)
}
