\name{pdnn.params.chiptype}
\alias{pdnn.params.chiptype}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A function to fit PDNN parameters }
\description{
  A function to fit PDNN parameters that are chip-type specific
}
\usage{
pdnn.params.chiptype(energy.param.file, probes.file = NULL, probes.pack= NULL,
                     probes.data.frame = NULL,
                     seq.name, x.name, y.name, affyid.name, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{energy.param.file}{ Path to the energy data file (see details) }
  \item{probes.file}{ Path to the probe files (see details) }
  \item{probes.pack}{ Name of the probe pack (see details) }
  \item{probes.data.frame}{ A \code{data.frame}}
  \item{seq.name, x.name, y.name, affyid.name}{ The names of the columns
    in the data.frame from \code{probes.pack} or \code{probes.file} for
  the probe sequences, the X positions, the Y positions and the probe
  set ID respectively}
  \item{verbose}{ verbosity (\code{TRUE} or \code{FALSE}) }
}
\details{
  The parameters \code{probes.file}, \code{probes.pack} and
  \code{probes.data.frame} are mutually exclusive.
  The function fits PDNN parameters that are specific to chip-types
  (hence specific to the probe sequences).
  It requires data files like the one found on Li Zhang's web page:
  (http://odin.mdacc.tmc.edu/~zhangli/PerfectMatch/)
  This should be computed once for all for a given chip type.
  Computed values for the chips are included in the package.
  This allows 'automagic' use of them when these chips types are used
  (as done in the function \code{\link{expressopdnn}}).
}
\value{
  A list of:
  \item{Eg}{environment. One entry per dinucleotide.}
  \item{Wg}{numerical vector}
  \item{En}{environment. One entry per dinucleotide.}
  \item{Wn}{numerical vector}
  \item{params.gene}{environment. One entry per gene, each entry is
    is a list of elements Sg, Sn, xy  and gene.i}
}
\section{warning}{ The X and Y positions in the \code{data.frame} are expected to be
  original ones in the Affymetrix files (starting at zero. They are offset by one within
  this function.}
\seealso{ \code{\link{find.params.pdnn}} }
\examples{
if (interactive()) {
  energy.file <- system.file("exampleData", "pdnn-energy-parameter_hg-u95av2.txt", package="affypdnn")

  params.chiptype <- pdnn.params.chiptype(energy.file, probes.pack="hgu95av2probe")
}

}
\keyword{ manip }

