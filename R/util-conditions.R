stop_no_method <- function(func, x, msg = NULL, msg_append = NULL) {
  if (is.null(msg))
    msg <- function(cls) paste0("no suitable ", substitute(func),
                                "() method for object of classes <",
                                paste0(cls, collapse = ", "), ">")
  msg <- msg(class(x))
  if (!is.null(msg_append))
    msg <- paste0(msg, "; ", msg_append)
  
  err <- structure(
    list(
      message = msg,
      call = NULL,
      func = substitute(func),
      classes = class(x)
    ),
    class = c("error_no_method", "error", "condition")
  )
  stop(err)
}
