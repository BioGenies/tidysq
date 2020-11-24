handle_warning_message <- function(msg,
                                   on_warning = getOption("tidysq_on_warning")) {
  if (!is.null(msg) && msg != "")
    switch(on_warning,
           error = stop(msg, call. = FALSE),
           warning = warning(msg, call. = FALSE),
           message = message(msg),
           none = invisible())
}
