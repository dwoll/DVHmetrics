\name{getEUD}
\alias{getEUD}
\alias{getEUD.DVHs}
\alias{getEUD.DVHLst}
\alias{getEUD.DVHLstLst}
\title{Generalized equivalent uniform dose (gEUD)}
\description{
Calculate generalized equivalent uniform dose (gEUD). May be based on EQD2.
}
\usage{
getEUD(x, EUDa, EUDfn=NULL, EUDab=NULL, ...)

\method{getEUD}{DVHs}(x, EUDa, EUDfn=NULL, EUDab=NULL, ...)

\method{getEUD}{DVHLst}(x, EUDa, EUDfn=NULL, EUDab=NULL, ...)

\method{getEUD}{DVHLstLst}(x, EUDa, EUDfn=NULL, EUDab=NULL, ...)
}
\arguments{
  \item{x}{One cumulative DVH (object of class \code{DVHs}, multiple cumulative DVHs from one patient with multiple structures (object of class \code{DVHLst}), or multiple cumulative DVHs from many patients, each with multiple structures (object of class \code{DVHLstLst}). See \code{\link[DVHmetrics]{readDVH}}.}
  \item{EUDa}{Exponential parameter a.}
  \item{EUDfn}{If gEUD should be based on EQD2: Number of fractions.}
  \item{EUDab}{If gEUD should be based on EQD2: alpha/beta ratio for the relevant tissue.}
  \item{...}{Ignored. Used to catch additional arguments passed from \code{\link[DVHmetrics]{getMetric}}.}
}
\value{
A data frame with variables \code{EUD}, \code{patID}, and \code{structure}.
}
\references{
Niemierko, A. (1999). A generalized concept of equivalent uniform dose. Medical Physics, 26(6), 1100.

Wu et al. (2002). Optimization of intensity modulated radiotherapy plans based on the equivalent uniform dose. International Journal of Radiation Oncology Biology Physics, 52, 224-235.
}
\seealso{
\code{\link[DVHmetrics]{getEQD2}},
\code{\link[DVHmetrics]{getMetric}}
}
\examples{
getEUD(dataMZ[[1]], EUDa=2)

# based on EQD2
getEUD(dataMZ[[1]], EUDa=2, EUDfn=20, EUDab=4)
}
