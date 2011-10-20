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
