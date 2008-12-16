expressopdnn <-  function(abatch,
                          ## --
                          bg.correct=FALSE,
                          bgcorrect.method = NULL,
                          bgcorrect.param = list(),
                          ## --
                          normalize = FALSE,
                          normalize.method = NULL,
                          normalize.param=list(),
                          ## --
                          pmcorrect.method = c("pdnn", "pdnnpredict"),
                          ##pmcorrect.param = list(),
                          ## --
                          findparams.param = list(),
                          summary.subset = NULL,
                          ## ---
                          eset.normalize = TRUE,
                          scale.to = 500,
                          verbose = TRUE
                     ) {

  nchips <- length(abatch)
  
###background stuff must be added before normalization!
  
  ## -- background correction method
  if (bg.correct && is.null(bgcorrect.method)) {
    stop("bg.method missing")
  }

  ## -- normalize.method
  if ((normalize) & (is.null(normalize.method))) {
    stop(paste("normalization method missing. Try one of:",
               normalize.methods(abatch), sep="\n"))
  }

  ## -- pm correction method
  if (is.null(pmcorrect.method)) {
    stop("pmcorrect.method missing")
  } else {
    pmcorrect.method <- match.arg(pmcorrect.method)
  }
  
  ## -- summary of what will be done
  if (verbose) {
    if (bg.correct){
      cat("background correction:", bgcorrect.method, "\n")
    }
    if (normalize) {
      cat("normalization:", normalize.method, "\n")
    }
    cat("PM/MM correction :", "pdnn", "\n")
    cat("expression values:", "pdnn", "\n")
  }
  
    ## -- background correct (if needed)
  if (bg.correct) {
    
    if (verbose)
      cat("background correcting...")
    
    abatch <- do.call(bg.correct, c(alist(abatch, method=bgcorrect.method),
                                       bgcorrect.param))
    
    if (verbose)
      cat("done.\n")
  }


  ## -- normalize (if wished)
  if (normalize) {
                                                                                
    if (verbose)
      cat("normalizing...")
    
    abatch <- do.call(normalize,
                      c(alist(abatch, normalize.method), normalize.param))
    
    if (verbose)
      cat("done.\n")
  }

  ## chip-type specific parameters
  params.chiptype <- findparams.param$params.chiptype
  if (is.null(params.chiptype)) {
    ## try to get it from the pack
    namebase <- cleancdfname(abatch@cdfName) 
    dataname <- paste(substr(namebase, 1,  nchar(namebase) - 3), ".pdnn.params", sep="")
    if(! dataname %in% do.call(data, list(package="affypdnn"))$results[, 3])
      stop("params.chiptype missing !")
    do.call(data, list(dataname, package="affypdnn"))
    assign("params.chiptype", get(dataname))
  }
  findparams.param$params.chiptype <- params.chiptype
  params <- do.call(find.params.pdnn,
                    c(alist(abatch), findparams.param))
  
  eset <- computeExprSet(abatch, pmcorrect.method=pmcorrect.method,
                         summary.method="pdnn",
                         ids=summary.subset,
                         pmcorrect.param = list(params,
                           params.chiptype=params.chiptype,
                           callingFromExpresso=TRUE),
                         summary.param = list(params))

  if (eset.normalize)
    eset <- pdnn.scalevalue.exprSet(eset, scale.to=scale.to)
  
  return(eset)
  
}
