.First.lib <- function(libname, pkgname, where) {
  
  require(affy)

  ## register the new summary method
  where <- match(paste("package:", pkgname, sep=""), search())
  cat("registering new summary method 'pdnn'.\n")
  assign("express.summary.stat.methods", c(express.summary.stat.methods, "pdnn"), envir=as.environment("package:affy"))
  cat("registering new pmcorrect method 'pdnn' and 'pdnnpredict'.\n")
  assign("pmcorrect.methods", c(pmcorrect.methods, "pdnn", "pdnnpredict"), envir=as.environment("package:affy"))
  ##assign("express.summary.stat.methods", c(express.summary.stat.methods, "pdnn"), envir=where)
  
}
