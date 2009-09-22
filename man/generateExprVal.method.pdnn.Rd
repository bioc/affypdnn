\name{generateExprVal.method.pdnn}
\alias{generateExprVal.method.pdnn}
\alias{pmcorrect.pdnn}
\alias{pmcorrect.pdnnpredict}
\title{Compute PM correction and summary expression value}
\description{
  Computes PM correction and summary expression value with PDNN method.
}
\usage{
pmcorrect.pdnn(object, params, gene=NULL, gene.i=NULL,
               params.chiptype=NULL, outlierlim=3, callingFromExpresso=FALSE)
pmcorrect.pdnnpredict(object, params, gene=NULL, gene.i=NULL,
               params.chiptype=NULL, outlierlim=3, callingFromExpresso=FALSE)
generateExprVal.method.pdnn(probes, params)
}
\arguments{
  \item{object}{object of \code{\link[affy:ProbeSet-class]{ProbeSet}}.}
  \item{probes}{matrix of PM-corrected signals (should be coming out of
  \code{pmcorrect.pdnn}).}
  \item{params}{experiments specific parameters.}
  \item{gene}{gene (probe set) ID (from wich the \code{gene.i} would be
  derived).}
  \item{gene.i}{gene index (see details).}
  \item{params.chiptype}{chip-specific parameters.}
  \item{outlierlim}{threshold for tagging a probe as an outlier.}
  \item{callingFromExpresso}{ is the function called through
    expresso. DO NOT play with that.}
}
\details{
  Only one of \code{gene}, \code{gene.i} should be specified. For most
  the users, this is \code{gene}.
  \code{\link{pmcorrect.pdnn}} and \code{\link{pmcorrect.pdnnpredict}}
  return what is called GSB and GSB + NSB + B in the paper by Zhang Li
  and collaborators.
}
\value{
  \code{pmcorrect.pdnn} and \code{pmcorrect.pdnnpredict} return a matrix (one row per probe, one column
  per chip) with attributes attached. \code{generateExprVal} returns a
  \code{list}:
  \item{exprs}{expression values}
  \item{se.exprs}{se expr. val.}
}
\seealso{\code{\link{pdnn.params.chiptype}}}
\examples{
data(hgu95av2.pdnn.params)
library(affydata)
data(Dilution)

## only one CEL to go faster
abatch <- Dilution[, 1]

## get the chip specific parameters
params <- find.params.pdnn(abatch, hgu95av2.pdnn.params)

## The thrill part: do we get like in the Figure 1-a of the reference ?
par(mfrow=c(2,2))
##ppset.name <- sample(featureNames(abatch), 2)
ppset.name <- c("41206_r_at", "31620_at")
ppset <- probeset(abatch, ppset.name)
for (i in 1:2) {
  ##ppset[[i]] <- transform(ppset[[i]], fun=log) # take the log as they do
  probes.pdnn <- pmcorrect.pdnnpredict(ppset[[i]], params,
                                       params.chiptype=hgu95av2.pdnn.params)
  ##probes.pdnn <- log(probes.pdnn)
  plot(ppset[[i]], main=paste(ppset.name[i], "\n(raw intensities)"))
  matplotProbesPDNN(probes.pdnn, main=paste(ppset.name[i], "\n(predicted intensities)"))
}

## pick the 50 first probeset IDs
## (to go faster)
ids <- featureNames(abatch)[1:100]

## compute the expression set (object of class 'ExpressionSet')
eset <- computeExprSet(abatch, pmcorrect.method="pdnn",
                       summary.method="pdnn", ids=ids,
                       summary.param = list(params, params.chiptype=hgu95av2.pdnn.params))

}
\keyword{manip}
