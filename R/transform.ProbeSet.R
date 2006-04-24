transform.ProbeSet <- function(`_data`, fun=I, ...) {
  `_data`@pm <- fun(`_data`@pm, ...)
  `_data`@mm <- fun(`_data`@mm, ...)
  return(`_data`)
}
