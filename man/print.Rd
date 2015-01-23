\name{print.DVHs}
\alias{print.DVHs}
\alias{print.DVHLst}
\alias{print.DVHLstLst}
\title{Print basic information about one or more DVHs}
\description{Print basic information (patients, structures, dose range) about one or more DVHs.}
\usage{
\method{print}{DVHs}(x, ...)

\method{print}{DVHLst}(x, ...)

\method{print}{DVHLstLst}(x, ...)
}
\arguments{
  \item{x}{A single DVH (object of class \code{DVHs}), multiple DVHs from one patient/structure (object of class \code{DVHLst}), or multiple DVHs from many patients/structures (object of class \code{DVHLstLst}). See \code{\link[DVHmetrics]{readDVH}}.}
  \item{...}{Further arguments: \code{print.DVHLst(x, verbose=TRUE)} prints more information about each DVH.}
}
\value{
Prints summary information about the DVHs.
}
\seealso{
\code{\link[DVHmetrics]{readDVH}}
}
\examples{
data(dataMZ)
print(dataMZ)
print(dataMZ, verbose=TRUE)
}
