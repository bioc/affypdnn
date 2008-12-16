find.params.pdnn <- function(abatch, params.chiptype=NULL, optim.method="BFGS", verbose=TRUE, give.warnings=TRUE) {
  ## chip-type specific parameters
  if (is.null(params.chiptype)) {
    ## try to get it from the pack
    namebase <- cleancdfname(abatch@cdfName) 
    dataname <- paste(substr(namebase, 1,  nchar(namebase) - 3), ".pdnn.params", sep="")
    if(! dataname %in% do.call(data, list(package="affypdnn"))$results[, 3])
      stop("params.chiptype missing !")
    do.call(data, list(dataname, package="affypdnn"))
    assign("params.chiptype", get(dataname))
  }
  if (verbose)
    cat("initializing data structure...")

  ## the following is a bit hard to follow (at least for the author of the
  ## code 8 months later ;) )... comments were added...

  ## names of probe sets (sorted)
  names.abatch <- sort(featureNames(abatch))
  ## number of probes in probes sets (set to zero at start)
  n.probes <- rep(0, length(names.abatch))
  ## number of CEL files in the abatch
  n.cel <- length(sampleNames(abatch))
  ##names.pdnnparams <- sort(ls(params.chiptype$params.gene))
  ## index for the names
  names.i <- rep(as.integer(NA), length(names.abatch))
  ## list of parameters 'lambda' (one list element for each probe set)  
  lambda <- vector("list", length=length(names.abatch))
  ## vector of parameters 'Fs' (one for each chip)
  Fs <- rep(as.numeric(NA), length=n.cel)
  ## vector of parameters 'S.lambda.E' (one for each probe set)
  S.lambda.E <- vector("numeric", length=length(names.abatch)) ## \sum \frac{\lambda_{ij}} {1+exp(E_{ij})}
  
  ##K1top <- vector("numeric", length=length(names.abatch))
  ##K2top <- vector("numeric", length=length(names.abatch))
  ##K3top <- vector("numeric", length=length(names.abatch))
  ##prealpha <- alpha <- beta <- rep(as.numeric(NA), length(names.abatch))
  Bs <- vector("numeric", length=n.cel)
  Ns <- vector("numeric", length=n.cel)

  ## init data structures
  for (gene.i in seq(along=names.abatch)) {
    gene <- names.abatch[gene.i]
    if(! exists(gene, envir=params.chiptype$params.gene)) {
      if (give.warnings) {
        warning(paste("The probeset", gene, "could not be found in the parameters (possible entanglement with cdfenvs, or missing probe sequence) !"))
      }
      n.probes[gene.i] <- 0
      names.abatch[gene.i] <- NA
    } else {
      names.i[gene.i] <- get(gene, envir=params.chiptype$params.gene)
      n.probes[gene.i] <- length(params.chiptype$gene.Sn[[names.i[gene.i]]])
      lambda[[gene.i]] <- matrix(as.numeric(NA), n.probes[gene.i], n.cel)
    }
    S.lambda.E[gene.i] <- NA
  }

  names.abatch.nonNA <- which(! is.na(names.abatch))
  i.pm <- indexProbes(abatch, which="pm", genenames=names.abatch)
  all.i.pm <- unlist(i.pm[names.abatch.nonNA])

  if (verbose)
    cat("done.\n")
  
  ## loop across chips
  for (cel.i in seq(along=sampleNames(abatch))) {
    
    if (verbose)
      cat("dealing with CEL", cel.i, ":\n")

    if (verbose)
      cat("  step 1...")
    
    for (gene.i in names.abatch.nonNA) {

      gene <- names.abatch[gene.i]
      
      gene.params.i <- names.i[gene.i]
      gene.Sg <- params.chiptype$gene.Sg[[gene.params.i]]
      gene.Sn <- params.chiptype$gene.Sn[[gene.params.i]]
      gene.xy <- params.chiptype$gene.xy[[gene.params.i]]
      gene.o <- match(i.pm[[gene.i]], xy2indices(gene.xy[,1], gene.xy[,2], abatch=abatch))
      probe.intensities <- intensity(abatch)[i.pm[[gene.i]], cel.i, drop=FALSE]
      
      lambda[[gene.i]][, cel.i] <- sqrt(probe.intensities * gene.Sg[gene.o])
      S.lambda.E[gene.i] <- sum( sqrt(probe.intensities / gene.Sg[gene.o]) )
      ##this.lambda <- lambda[[gene.i]][, cel.i]
      ##lambda[[gene]] <- sqrt(intensity(abatch)[i.pm, ] * get(gene, envir=params.chiptype$Sg))
       ## prealpha[gene.i] <- sum(1 / (gene.Sg[gene.o] * this.lambda))
      ##K1top[gene.i] <- sum(probe.intensities / this.lambda)
      ##K2top[gene.i] <- sum(1 / this.lambda)
      ##K3top[gene.i] <- sum(1/(gene.Sn[gene.o] * this.lambda))
    }
    
    if (verbose)
      cat("done.\n")

    ##K1top.rep <- rep(K1top, n.probes)
    ##K2top.rep <- rep(K2top, n.probes)
    ##K3top.rep <- rep(K3top, n.probes)
    
    ##all.gene.params <- multiget(names.abatch, envir=params.chiptype$params.gene)
    ##alpha <- rep(prealpha, n.probes) * unlist(lapply(all.gene.params, function(x) x$Sg))
    ##KI1 <- (K1top.rep / alpha) / intensity(abatch)[all.i.pm, cel.i, drop=FALSE]
    ##KI2 <- (K2top.rep / alpha - 1) / intensity(abatch)[all.i.pm, cel.i, drop=FALSE]
    ##KI3 <- (K3top.rep / alpha - (1 / unlist(lapply(all.gene.params, function(x) x$Sn)))) /
    ##  intensity(abatch)[all.i.pm, cel.i, drop=FALSE]
    ##rm(all.gene.params, K1top.rep, K2top.rep, K3top.rep)
    
    if (verbose)
      cat("  step 2...")
    ##FIXME: starting values always the same ?
    B<-50
    N<-3000

    lambda.cel <- unlist(lapply(names.abatch.nonNA, function(x) lambda[[x]][,cel.i]))
    Betaij.left <- sum(1 / unlist(params.chiptype$gene.Sg[seq(along=names.abatch.nonNA)]) / lambda.cel)
    Betaij <- lapply(seq(along=names.abatch.nonNA), function(i, x, y) y[[i]] * x, Betaij.left, params.chiptype$gene.Sg)
    Gij.top <- sum(intensity(abatch)[all.i.pm, cel.i] / lambda.cel)
    Gij <- lapply(seq(along=names.abatch.nonNA), function(i, x, y) x / y[[i]], Gij.top, Betaij)
    Hij.top <- sum(1 / lambda.cel)
    Hij <- lapply(seq(along=names.abatch.nonNA), function(i, x, y) x / y[[i]] - 1, Hij.top, Betaij)
    Kij.lefttop <- sum(1 / lambda.cel / unlist(params.chiptype$gene.Sn[seq(along=names.abatch.nonNA)]))
    Kij <- lapply(seq(along=names.abatch.nonNA), function(i, x, y, z) x / y[[i]] + 1 / z[[i]], Kij.lefttop, Betaij, params.chiptype$gene.Sn)
    intensity.forfit <- intensity(abatch)[unlist(i.pm[names.abatch.nonNA]), cel.i]
    intensity.forfit.positive <- intensity.forfit > 0
    Gij.forfit <- unlist(Gij[names.i[names.abatch.nonNA]])
    Hij.forfit <- unlist(Hij[names.i[names.abatch.nonNA]])
    Kij.forfit <- unlist(Kij[names.i[names.abatch.nonNA]])
    
    ##fit.f<-function(par, Gij, Hij, Kij, intensities.matrix, i.pm){
    ##fit.f<-function(par){
    fit.f<-function(par, Gij.forfit, Hij.forfit, Kij.forfit) {
      B <- par[1]
      N <- par[2]
      sum.sqrt <- 0
      n.probes <- 0
      Ihat <-  Gij.forfit - B * Hij.forfit - N * Kij.forfit
      good <- Ihat > 0 & intensity.forfit.positive
      sum.sqr <- sum((log(Ihat[good] / intensity.forfit[good]))^2)
#       for (gene.i in names.abatch.nonNA) {
#         ##cat(gene.i, "\n")
#         i <- names.i[gene.i]
#         probe.intensities <- intensity.matrix[i.pm[[gene.i]], cel.i]
#         I.hat <- Gij[[i]] - B * Hij[[i]] - N * Kij[[i]]
#         good <- I.hat > 0 & probe.intensities > 0
#         n.probes <- n.probes + sum(good)
#         sum.sqr <- sum.sqrt + sum((log(I.hat[good]) - log(probe.intensities[good]))^2)
#       }
      ##cat(sum.sqr / sum(good), "\n")
      return(sum.sqr / sum(good))
      ##return( sum(log(Y[Y >= 0])**2) / length(Y[Y >= 0]))
    }

    if (optim.method == FALSE) {
      ## use steepest descent
      i<-succes<-failed<-0

      DeltaN<-500;DeltaB<-10
      Fit <- +Inf
      Bdirection<-1
      Ndirection<-1
      improvement<-2
      
      while (improvement>1.00000001 || DeltaN>0.5 || DeltaB>0.5){
        DB<-DeltaB*runif(1,0.85,1)
        Bnew<-B+Bdirection*DB
        DN<-DeltaN*runif(1,0.85,1)
        Nnew<-N+Ndirection*DN
        Fnew <- fit.f(c(Bnew,Nnew), Gij.forfit=Gij.forfit, Hij.forfit=Hij.forfit, Kij.forfit=Kij.forfit)
        i<-i+1
        if(Fnew<Fit){# Succes
                                        #cat("PDNN: Iteration:",i, " N':", N, " B:", B, " F:", Fnew, "\n")#(Ndirection*DeltaN)/(Bdirection*DeltaB)
          improvement <- Fit/Fnew
          succes<-succes+1
          failed<-0
          DeltaB<-DB; DeltaN<-DN
          Fit<-Fnew;
          B<-Bnew;
          N<-Nnew
          if(succes>1){
            increase<-runif(1,1,1.5)
            DeltaB<-DeltaB*increase
            DeltaN<-DeltaN*increase
          }
          
        }else{ # bad guess
          failed<-failed+1
          succes<-0
          if(failed>1){
            DeltaB<-DeltaB*runif(1,0.85**failed,1)
            DeltaN<-DeltaN*runif(1,0.9**failed,1)
            if(runif(1,0,1)>0.4999){Bdirection<-1}else{Bdirection<--1}
            if(runif(1,0,1)>0.4999){Ndirection<-1}else{Ndirection<--1}
          }
        }
      }
      Bs[cel.i] <- B
      Ns[cel.i] <- N
    } else {
      ## use R's 'optim'
      
      ##par.f <- optim(c(B, N), fit.f, method=optim.method, lower=-Inf, upper=-Inf, control=list(), hessian=FALSE, Betaij, Gij, Hij, Kij)
      par.f <- optim(c(B, N), fit.f, method=optim.method, Gij=Gij.forfit, Hij=Hij.forfit, Kij=Kij.forfit)
      
      Bs[cel.i] <- par.f$par[1]
      Ns[cel.i] <- par.f$par[2]
    }
    if (verbose)
      cat("done.\n")

    Fs[cel.i] <- fit.f(c(Bs[cel.i], Ns[cel.i]), Gij.forfit, Hij.forfit, Kij.forfit)
    
  }
  
  return(list(lambda=lambda, Bs=Bs, Ns=Ns, Fs=Fs, names.abatch=names.abatch, names.i=names.i))
}
