\name{transform.ProbeSet}
\alias{transform.ProbeSet}
\title{ A function to tranform a ProbeSet }
\description{
  A function to transform the PMs and MMs in a ProbeSet
}
\usage{
transform.ProbeSet(`_data`, fun = I, ...)
}
\arguments{
  \item{_data}{ object of \code{\link[affy]{ProbeSet-class}} }
  \item{fun}{ a function. The identity function by default. }
  \item{\dots}{ optional arguments for \code{fun}}
}
\details{
  The function \code{fun} is applied to the slots \code{pm} and
  \code{mm}.
  The function \code{vsnh} in the package \code{vsn} is a recommended
  argument for \code{fun}.
}
\value{
  An object of class \code{ProbeSet}.
}
\author{ Laurent }
\note{ This function should make to the package \code{affy} for the version
  1.4.x.}
\examples{
library(affydata)

data(Dilution)

ppset.name <- sample(geneNames(Dilution), 1)
ppset <- probeset(Dilution, ppset.name)[[1]]
ppset.log <- transform(ppset, fun=log)

par(mfrow=c(1,2))
plot(ppset)
plot(ppset.log)

}
\keyword{ manip }

