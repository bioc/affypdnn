\name{find.params.pdnn}
\alias{find.params.pdnn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A function to find the experiment specific PDNN parameters}
\description{
  A function to find the paramaters specific to the chips in an
  AffyBatch object.
}
\usage{
find.params.pdnn(abatch, params.chiptype, optim.method = "BFGS", verbose = TRUE, give.warnings=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{abatch}{ an instance of \code{\link[affy]{AffyBatch-class}}. }
  \item{params.chiptype}{ chip-type specific parameters (see details)}
  \item{optim.method}{method for the optimization function
    \code{\link{optim}}. When \code{FALSE}, a steepest-descent method of our own
  is used.}
  \item{verbose}{ verbosity (\code{TRUE} or \code{FALSE}) }
  \item{give.warnings}{ report probeset IDs in the \code{abatch} that could not be found in
    the \code{params.chiptype}}
}
\details{
  This function fits PDNN parameters that are specific to experimental
  values. The parameters common to all the chips of a certain type
  are returned by the function \code{\link{pdnn.params.chiptype}}.
  If \code{NULL}, the parameter files included in the
  package will be used whenever possible...
}
\value{
  A list of
  \item{lambda}{The lambda's}
  \item{Bs}{The B's}
  \item{Ns}{The N's}
  \item{Fs}{The F's}
}
\references{ Li Zhang, Michael F. Miles and Kenneth D. Aldape - A model of molecular interactions on short oligonucleotide
  arrays, 2003, Nature Biotech., vol. 21, n.7 }

\seealso{ \code{\link{pdnn.params.chiptype}}, \code{\link{generateExprVal.method.pdnn}}}
\examples{
## load a chip-specific parameter file
## (as returned by the function pdnn.params.chiptype)
data(hgu95av2.pdnn.params)

## load experimental data
library(affydata)
data(Dilution)

## one CEL to go faster
afbatch <- Dilution[, 1]
params <- find.params.pdnn(afbatch, hgu95av2.pdnn.params, optim.method =
FALSE, give.warnings=FALSE)

}
\keyword{ manip }

