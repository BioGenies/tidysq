stop_no_method <- function(func, x) {
  err <- structure(
    list(
      message = paste0("no suitable ", substitute(func),
                       "() method for object of classes <",
                       paste0(class(x), collapse = ", "), ">"),
      call = NULL,
      func = substitute(func),
      classes = class(x)
    ),
    class = c("error_no_method", "error", "condition")
  )
  stop(err)
}
