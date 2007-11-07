matplotProbesPDNN <- function(x, type="l", ...) {
  matplot(x, type=type, ...)  
  ok <- attr(x, "ok")
  for (i in seq(1, ncol(x), length=ncol(x))) {
    points(x[, i], pch=c(1,3)[as.integer(ok[, i])+1])
  }
}

pmcorrect.pdnnpredict <- function(object, params, gene=NULL, gene.i=NULL, params.chiptype=NULL, outlierlim=3, callingFromExpresso=FALSE) {  

  new.int <- pmcorrect.pdnn(object, params, gene=gene, gene.i=gene.i, params.chiptype=params.chiptype, outlierlim=outlierlim, callingFromExpresso=callingFromExpresso)
  ## do the prediction bit
  
  Bs <- params$Bs
  Ns <- params$Ns
  Sn.gene <- attr(new.int, "Sn.gene")

 
  new.int <- sweep(new.int + 1/outer(Sn.gene, Ns, "/"), 2, Bs, "+")
  
  return(new.int)
}

pmcorrect.pdnn <- function(object, params, gene=NULL, gene.i=NULL, params.chiptype=NULL, outlierlim=3, callingFromExpresso=FALSE) {  

  probes <- object@pm
  if (is.null(params.chiptype)) {
    stop("params.chiptype must be specified.")
  }

  lambda <- params$lambda
  Bs <- params$Bs
  Ns <- params$Ns
  Fs <- params$Fs

  if (is.null(gene.i)) {
    if (callingFromExpresso) {
      gene <- get("id", envir=parent.frame(3)) ## dynamic lexical scoping... (not static)
    } else {
      if (is.null(gene))
        ##stop("gene.i must be specified.")
        gene <- object@id
    }
  }

  gene.i <- get(gene, envir=params.chiptype$params.gene)
  Sg.gene <- params.chiptype$gene.Sg[[gene.i]]
  Sn.gene <- params.chiptype$gene.Sn[[gene.i]]
  
  Ntop <- sweep(probes, 2, Bs , "-") - sapply(Ns, "/", Sn.gene)
  
  new.int <- probes
  ok <- matrix(as.logical(NA), nr=nrow(probes), nc=ncol(probes))
  
  for (cel.i in seq(1, ncol(probes), length=ncol(probes))) {
    new.int[, cel.i] <- sum((Ntop / lambda[[gene.i]][, cel.i]) / sum(1 / Sg.gene / lambda[[gene.i]][, cel.i])) /
      Sg.gene #+ Ns[cel.i] / Sn.gene + Bs[cel.i]
    
    ok[, cel.i] <- Ntop[, cel.i] >= 0 & probes[, cel.i] / new.int[, cel.i] > 0 & log(probes[, cel.i] / new.int[, cel.i]) < outlierlim * sqrt(Fs[cel.i]) & !is.na(log(probes[, cel.i] / new.int[, cel.i]))
  }
  attr(new.int, "ok") <- ok
  attr(new.int, "Ntop") <- Ntop
  attr(new.int, "gene.i") <- gene.i
  attr(new.int, "Sg.gene") <- Sg.gene
  attr(new.int, "Sn.gene") <- Sn.gene

  return(new.int)
}


generateExprVal.method.pdnn <- function(probes, params) {  
  
#   if (is.null(params.chiptype)) {
#     stop("params.chiptype must be specified.")
#   }

  if (is.null(attr(probes, "ok")) || is.null(attr(probes, "Ntop")) || is.null(attr(probes, "gene.i")))
    stop("The summary method 'pdnn' can only be used with the pmcorrect method 'pdnn' !")

  lambda <- params$lambda
#   Bs <- params$Bs
#   Ns <- params$Ns
#   Fs <- params$Fs

  
#   if (is.null(gene.i)) {
#     if (callingFromExpresso)
#       gene <- get("id", envir=parent.frame(3)) ## dynamic lexical scoping... (not static)
#     else
#       stop("gene.i must be specified.")
#   }    
#   params.gene <- get(gene, envir=params.chiptype$params.gene)
#   gene.i <- params.gene$gene.i
#   Sg.gene <- params.gene$Sg
#   Sn.gene <- params.gene$Sn

  ##gene <- names.abatch[gene.i]
  ##names(Nj) <- names(Iobs)
  ##options(warn=-1)

  ##FIXME: useless ?
  ##i.pm <- indexProbes(abatch, gene, which="pm")
  ##i.mm <- indexProbes(abatch, gene, which="mm")

  ##cat(str(probes))
  ##cat(str(Bs+Ns))
  ##cat(str(get(gene, envir=params.chiptype$Sn)))
  
#   Ntop <- sweep(probes, 2, Bs , "-") - sapply(Ns, "/", params.gene$Sn)
  
#   new.int <- probes
#   ok <- matrix(as.logical(NA), nr=nrow(probes), nc=ncol(probes))
  
#   for (cel.i in seq(1, ncol(probes), length=ncol(probes))) {
#     new.int[, cel.i] <- sum((Ntop / lambda[[gene.i]][, cel.i]) / sum(1 / Sg.gene / lambda[[gene.i]][, cel.i])) /
#       Sg.gene + Ns[cel.i] / Sn.gene + Bs[cel.i]
#     ok[, cel.i] <- Ntop[, cel.i] >= 0 & probes[, cel.i] / new.int[, cel.i] > 0 & log(probes[, cel.i] / new.int[, cel.i]) < outlierlim * sqrt(Fs[cel.i]) & !is.na(log(probes[, cel.i] / new.int[, cel.i]))
#   }

  Ntop <- attr(probes, "Ntop")
  Sg.gene <- attr(probes, "Sg.gene")
  gene.i <- attr(probes, "gene.i")
  ok <- attr(probes, "ok")
  
  expr.val <- rep(as.numeric(NA), ncol(probes))
  expr.se <- rep(as.numeric(NA), ncol(probes))
  ##cat(gene.i, "-- ")
  for (cel.i in seq(1, ncol(probes), length=ncol(probes))) {
    sum((Ntop / lambda[[gene.i]][, cel.i])[ok[, cel.i]]) / sum((1/Sg.gene / lambda[[gene.i]][, cel.i])[ok[, cel.i]])
    expr.val[cel.i] <- sum((Ntop / lambda[[gene.i]][, cel.i])[ok[, cel.i]]) / sum((1/Sg.gene / lambda[[gene.i]][, cel.i])[ok[, cel.i]])
    
  }
  return(list(exprs=expr.val, se.exprs=expr.se))
}

pdnn.scalevalue.exprSet <- function(eset, scale.to=500) {
  m <- exprs(eset)
  m.mean <- apply(m, 2, mean, na.rm=TRUE)
  mm <- sweep(m, 2, scale.to/m.mean, "*")
  exprs(eset) <- mm
  return(eset)
}
