\name{hgu133a.pdnn.params}
\alias{chiptype.pdnn.params}
%\alias{hgu133a.pdnn.params}
\alias{hgu95av2.pdnn.params}
%\alias{mgu74av2.pdnn.params}
%\alias{pdnn-energy-parameter_hg_u133a}
\alias{pdnn-energy-parameter_hg-u95av2}
%\alias{pdnn-energy-parameter_mg_u74av2}
\docType{data}
\title{Chip-type specific data}
\description{
Chip-type specific data structure.
}
\usage{
%data(hgu133a.pdnn.params)
data(hgu95av2.pdnn.params)
}
\format{
  The format is a list:
\describe{
  \item{Eg}{environment (one entry per dinucleotide)}
  \item{Wg}{numerical vector.}
  \item{En}{environment (one entry per dinucleotide)}
  \item{Wn}{numerical vector.}
  \item{gene.Sn}{list (one entry per affyID)}
  \item{gene.Sg}{list (one entry per affyID)}
  \item{gene.xy}{list (one entry per affyID)}
  \item{params.gene}{environment}
}
}
\details{
  These chip-specific data structures are generated from the data files
  made available by the author of the PDNN paper (see the section source). They
  are stored as \code{data} to save some computation time.
  The data structures were made using the function
  \code{\link{pdnn.params.chiptype}}. The data files are included in the
  the directory `data' of the package. 
}
\note{
  To lower the size of the package, the only chip-specific data structures
  included in the package is the one for HG-U95Av2.
}
\source{
    Li Zhang, Michael F. Miles and Kenneth D. Aldape - A model of molecular interactions on short oligonucleotide
  arrays, 2003, Nature Biotech., vol. 21, n.7
}
\examples{
## give the path the original energy parameter files included in the package
list.files(system.file("exampleData", package="affypdnn"),
           "^pdnn-energy-parameter_", full.names=TRUE)
}
\keyword{datasets}
