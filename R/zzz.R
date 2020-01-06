.onLoad <- function(libname, pkgname) {
  
  require(affy)

  ## register the new summary method
  where <- match(paste("package:", pkgname, sep=""), search())
  cat("registering new summary method 'pdnn'.\n")
  upDate.express.summary.stat.methods(c(affy::express.summary.stat.methods(), 
        "pdnn"))
  cat("registering new pmcorrect method 'pdnn' and 'pdnnpredict'.\n")
  upDate.pmcorrect.methods(c(affy::pmcorrect.methods(), "pdnn", "pdnnpredict"))
}
.onAttach <- function(libname, pkgname) {
    msg <- sprintf(
        "Package '%s' is deprecated and will be removed from Bioconductor
         version %s", pkgname, "3.12")
    .Deprecated(msg=paste(strwrap(msg, exdent=2), collapse="\n"))
}
