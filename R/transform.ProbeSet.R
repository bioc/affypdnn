transform.ProbeSet <- function(x, fun=I, ...) {
  x@pm <- fun(x@pm, ...)
  x@mm <- fun(x@mm, ...)
  return(x)
}
