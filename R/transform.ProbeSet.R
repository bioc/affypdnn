transform.ProbeSet <- function(`_data`, fun=I, ...) {
  x <- `_data`
  x@pm <- fun(x@pm, ...)
  x@mm <- fun(x@mm, ...)
  return(x)
}
