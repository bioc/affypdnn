\name{matplotProbesPDNN}
\alias{matplotProbesPDNN}
\title{ Plot the PDNN computed probe intensities }
\description{
  Plot the probe intensities as computed by 'pmcorrect.pdnn' or 'pmcorrect.pdnnpredict'
}
\usage{
matplotProbesPDNN(x, type="l", ...)
}
\arguments{
  \item{x}{ a matrix (and attributes) as returned by \code{pmcorrect.pdnn}
    or \code{pmcorrect.pdnnpredict}. }
  \item{type}{type of plot (same as in \code{matplot})}
  \item{\dots}{ optional arguments to be passed to \code{matplot} }
}
\details{
  The crosses are the probe intensities which are considered `ok' by the
  outlier detection part of the algorithm, while the circles are the
  ones considered `outliers'
}
\value{
  Only used for its side-effect.
}
\seealso{\code{\link{pmcorrect.pdnn}} and \code{\link{pmcorrect.pdnnpredict}}}
\examples{
# see 'pmcorrect.pdnn'
}
\keyword{ hplot }
