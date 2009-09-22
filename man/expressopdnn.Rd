\name{expressopdnn}
\alias{affypdnn}
\alias{expressopdnn}
\alias{pdnn.scalevalue.exprSet}
\title{ Position Dependant Nearest Neighbors model for affy}
\description{
  A wrapper to perform the PDNN method.
}
\usage{
pdnn.scalevalue.exprSet(eset, scale.to=500)
expressopdnn(abatch,
      # background correction
             bg.correct = FALSE,
             bgcorrect.method = NULL,
             bgcorrect.param = list(),
      # normalize
             normalize = FALSE,
             normalize.method = NULL,
             normalize.param = list(),
%       # pm correction
              pmcorrect.method = c("pdnn", "pdnnpredict"),
%              pmcorrect.param = list(),
      # pdnn
             findparams.param = list(),          
      # expression values
             summary.subset = NULL,
      # PDNN expression values scaling
             eset.normalize = TRUE,
             scale.to = 500,
      # misc.
             verbose = TRUE)
}
\arguments{
  \item{abatch}{object of \code{\link[affy]{AffyBatch-class}}.}
  \item{bg.correct}{a boolean to express whether background correction
    is wanted or not.}
  \item{bgcorrect.method}{the name of the background adjustment method.}
  \item{bgcorrect.param}{a list of parameters for bgcorrect.method (if
    needed/wanted).}
  \item{eset}{an object of
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet-class}}.} 
  \item{normalize}{normalization step wished or not.}
  \item{normalize.method}{the normalization method to use.}
  \item{normalize.param}{a list of parameters to be passed to the
    normalization method (if wanted).}
   \item{pmcorrect.method}{the name of the PM adjustement method (only
     two choices here, default to 'pdnn').}
%   \item{pmcorrect.param}{a list of parameters for pmcorrect.method (if
%     needed/wanted)}
  \item{findparams.param}{a list of parameters to be passed to
    \code{find.params.pdnn}.} 
  \item{eset.normalize}{is any normalization step on expression values
    to be performed.}
  \item{scale.to}{a value to scale against.}
  \item{summary.subset}{a list of 'affyids'. If \code{NULL}, then an
    expression summary value is computed for everything on the chip.}
  \item{verbose}{logical value. If \code{TRUE} it writes out some
    messages.}
}
\details{
  \code{expressopdnn} is very similar to \code{\link[affy]{expresso}}. It is
  mainly a wrapper around the pre-processing steps `background
  correction', `normalization', `perfect match correction' and
  the PDNN method to compute expression values (see the first
  reference for more details about the preprocessing steps and
  and the second reference for further details about the PDNN method).

  The wrapper \code{expresso} has no way to handle easily the
  computation of chip-wide results that have to be used during
  the \code{computeExprSet} step. An easy way to overcome this
  was to write this simple wrapper.

  \code{pdnn.scalevalue} is performed after the expression values
  have computed to somehow `normalize' the values between different
  chips. When setting \code{normalize} to \code{TRUE} this step
  might be considered unnecessary (and the \code{eset.normalize} set
  to \code{FALSE}).
}
\value{
  An object of
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet-class}}, with
  an attribute 
  \code{pps.warnings} as returned by the method
  \code{\link[affy:generateExprSet-methods]{computeExprSet}}. 
}
\seealso{\code{\link[affy]{expresso}} and
  \code{\link{generateExprVal.method.pdnn}}} 
\examples{
## load pre-computed parameters
data(hgu95av2.pdnn.params)

library(affydata)
data(Dilution)

## one CEL to go faster
afbatch <- Dilution[, 1]

## Take only few IDs (the 10 first)
ids <- ls(getCdfInfo(afbatch))[1:10]
eset <- expressopdnn(afbatch, bg.correct=FALSE,
                     normalize=FALSE,
                     findparams.param=list(params.chiptype=hgu95av2.pdnn.params,
                                           give.warnings=FALSE),
                     summary.subset=ids)
}
\keyword{manip}
