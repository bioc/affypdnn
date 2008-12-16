pdnn.params.chiptype <- function(energy.param.file, probes.file = NULL, probes.pack = NULL, probes.data.frame = NULL,
                                 seq.name="sequence", x.name="x", y.name="y", affyid.name="Probe.Set.Name",
                                 verbose=TRUE) {


  if (! file.exists(energy.param.file))
    stop(paste("Cannot find the energy.param.file '", energy.param.file, sep=""))
  
  if (sum(c(is.null(probes.file), is.null(probes.pack), is.null(probes.data.frame))) != 2) {
    stop("Specify one 'probe.file' _or_ one 'probe.pack' _or_ 'probe.data.frame'")
  }

#   if (is.null(cdfName)) {
#     stop("'cdfName', a name to find the corresponding cdfenv is missing !")
#   }

#   ## Hack for version 1.2.x of 'affy'
#   a <- new("AffyBatch", cdfName=cdfName)
#   cdfenv <- getCdfInfo(a)
#   rm(a)
  
  ## FIXME
  if (!is.null(probes.file)) {
    probe.tab <- read.table(probes.file, sep="\t",header=TRUE, nrows=2)##[c(1,2,3,5)]
  }
  
  if (!is.null(probes.pack)) {
    do.call(library, list(probes.pack))
    probe.tab <- get(probes.pack, envir=as.environment(paste("package:", probes.pack, sep="")))
  }

  if (!is.null(probes.data.frame))
    probe.tab <- probes.data.frame
  ##

  xy.offset <- getOption("BioC")$affy$xy.offset
  
  i.seq <- match(seq.name, names(probe.tab))
  i.x <- match(x.name, names(probe.tab))
  i.y <- match(y.name, names(probe.tab))
  i.affyid <- match(affyid.name, names(probe.tab))

  if (any(is.na(c(i.seq, i.x, i.y, i.affyid)))) {
    m <- paste("You asked :", paste(" affyid.name:", affyid.name), paste(" x.name  :", x.name),
               paste(" y.name  :", y.name), paste(" seq.name :", seq.name),
               paste("While the names in the data.frame are:"),
               paste(paste("",names(probe.tab)), collapse="\n"), sep="\n")
    stop(m)
  }

  if (! is.null(probes.file)) {
    if (verbose)
      cat("Reading the probes file...")
    probe.tab <- read.table(probes.file, sep="\t",header=TRUE)
    if (verbose)
      cat("done.\n")
  }
  
  probe.x <- probe.tab[[i.x]] + xy.offset
  probe.y <- probe.tab[[i.y]] + xy.offset
  affy.id <- probe.tab[[i.affyid]]
  probe.seq <- tolower(as.character(probe.tab[[i.seq]]))
  
  if(verbose)
    cat("Calculating chip type specific parameters, (may take some time)...\n")

  ## Reading probe sequence, energy and position weight files
  
  ## FIXME (automagic to do) ?
  
  ep <- read.table(energy.param.file, nrows=80, header=TRUE, as.is=TRUE)

  Wg <- as.vector(ep[33:56, 2]) ## weights (specific)
  
  Wn <- as.vector(ep[57:80, 2]) ## weights (non-specific)

#    params.chiptype <- list(Eg = new.env(hash=TRUE),
#                            Wg = Wg,
#                            En= new.env(hash=TRUE),
#                            Wn = Wn,
#                            Sg = new.env(hash=TRUE),
#                            Sn = new.env(hash=TRUE),
#                            gene2i = new.env(hash=TRUE))
  
  params.chiptype <- list(Eg = new.env(hash=TRUE),
                          Wg = Wg,
                          En= new.env(hash=TRUE),
                          Wn = Wn,
                          gene.Sn = NA,
                          gene.Sg = NA,
                          gene.xy = NA,
                          gene.name = NA, ## only here for cross-checks
                          params.gene = new.env(hash=TRUE))
  
  
  Eg <- multiassign(as.list(ep[[1]][1:16]), as.list(as.vector(ep[1:16, 2])), envir=params.chiptype$Eg) ## energy (specific)
  
  ##Wg <- multiassign(rownames(ep), as.vector(ep[33:56, 2]), envir=env.list$Wg) ## weights (specific)
  
  En <- multiassign(as.list(ep[[1]][17:32]), as.list(as.vector(ep[17:32, 2])), envir=params.chiptype$En) ## energy (non-specific)
  
  
  ##Wg <- params.chiptype$Wg
  ##Wn <- params.chiptype$Wn
  En <- params.chiptype$En
  Eg <- params.chiptype$Eg
  
  Sf <- function(oligo){
    ## oligo: vector of probe sequences
    ## calculating (1/exp(E)) and (1/exp(E*)

    ENv <-EGv <-vector(length = length(oligo))

    increment <- as.integer(1)
    ## loop across the oligos
    for (g in 1:length(oligo)){
      
      EG <- EN <- 0

      ## walk along the sequence
      for (k in seq(1, nchar(oligo[g])-1)) {

        ## FIXME: build hashtables for Eg and En
        ##di.nucl <- substr(oligo[g], k, k+1)
        di.nucl <- .Internal(substr(oligo[g], k, k+increment))
        ## FIXME cast necessary ?
        ##EG <- as.numeric(EG + Wg[k] * get(di.nucl, envir = Eg))
        ##EN <- as.numeric(EN + Wn[k] * get(di.nucl, envir = En))
        EG <- EG + Wg[k] * get(di.nucl, envir = Eg)
        EN <- EN + Wn[k] * get(di.nucl, envir = En)
        
      }
      
      EGv[g] <- 1 + exp(EG)
      ENv[g] <- 1 + exp(EN)
    }
    
    return(cbind(EGv, ENv))
  }

  
  ## FIXME (speedup w/ affy cdfenvs ?)
  ##S <- tapply(as.vector(probe.tab[[i.seq]]), affy.id, Sf)
  S.index <- tapply(seq(1, nrow(probe.tab), length=nrow(probe.tab)), affy.id, I, simplify=FALSE)
  
  params.chiptype$gene.Sn <- vector("list", length=length(S.index))
  params.chiptype$gene.Sg <- vector("list", length=length(S.index))
  params.chiptype$gene.xy <- vector("list", length=length(S.index))
  params.chiptype$gene.name <- vector("character", length=length(S.index))
  
  params.gene <- params.chiptype$params.gene

  if (verbose) {
    pbt <- new("ProgressBarText", length(S.index), barsteps = as.integer(20))
    open(pbt)
  }
  
  ## FIXME:
  for (gene.i in seq(along=S.index)) {

    if (verbose)
      update(pbt)
    
    gene <- names(S.index)[gene.i]
    gene.S.index <- S.index[[gene.i]]
    sequences.in.ppset <- probe.seq[gene.S.index]
    S <- Sf(sequences.in.ppset)
    ## FIXME store the Xs and Ys for now (better scheme when cdfenvs in classes)
    
    ##assign(gene, S[[gene.i]][, 1], envir=Sg) #Sg is now (1+exp(E))
    ##assign(gene, S[[gene.i]][, 2], envir=Sn) #Sn is now (1+exp(E*))
    ##assign(gene, gene.i, envir=params.chiptype$gene2i)
    assign(gene, gene.i, envir=params.gene)
    params.chiptype$gene.Sg[[gene.i]] <- S[, 1]
    params.chiptype$gene.Sn[[gene.i]] <- S[, 2]
    params.chiptype$gene.xy[[gene.i]] <- cbind(probe.x[gene.S.index], probe.y[gene.S.index])
    params.chiptype$gene.name[[gene.i]] <- gene
  }
  
  if (verbose)
    close(pbt)
  
  return(params.chiptype)
  
}

