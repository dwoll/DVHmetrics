\name{getBED}
\alias{getBED}
\alias{getBED.default}
\alias{getBED.DVHs}
\alias{getBED.DVHLst}
\alias{getBED.DVHLstLst}
\title{Calculate biologically effective dose (BED)}
\description{
Calculate biologically effective dose (BED) according to the linear-quadratic model.
}
\usage{
getBED(D=NULL, fd=NULL, fn=NULL, ab=NULL)

\method{getBED}{default}(D=NULL, fd=NULL, fn=NULL, ab=NULL)

\method{getBED}{DVHs}(D=NULL, fd=NULL, fn=NULL, ab=NULL)

\method{getBED}{DVHLst}(D=NULL, fd=NULL, fn=NULL, ab=NULL)

\method{getBED}{DVHLstLst}(D=NULL, fd=NULL, fn=NULL, ab=NULL)
}
\arguments{
  \item{D}{Default: Total dose. If \code{NULL}, \code{fn} must be given. Alternative: One cumulative DVH (object of class \code{DVHs}, multiple cumulative DVHs from one patient with multiple structures (object of class \code{DVHLst}), or multiple cumulative DVHs from many patients, each with multiple structures (object of class \code{DVHLstLst}). See \code{\link[DVHmetrics]{readDVH}}. If some kind of \code{DVH} object, \code{fn} must be given.}
  \item{fd}{Fractional dose.}
  \item{fn}{Number of fractions. If \code{NULL}, \code{D} must be the total dose.}
  \item{ab}{alpha/beta ratio for the relevant tissue.}
}
\value{
If \code{D} is the total dose: a data frame with variables \code{BED} and \code{fractDose}.

If \code{D} is some kind of \code{DVH} object, a data frame with variables \code{BED}, \code{fractDose}, \code{patID}, and \code{structure}.
}
\references{
Fowler, J. F. (2010). 21 years of Biologically Effective Dose. British Journal of Radiology, 83, 554-568.
}
\seealso{
\code{\link[DVHmetrics]{getEQD2}},
\code{\link[DVHmetrics]{getIsoEffD}}
}
\examples{
getBED(D=50, fd=2.5, ab=c(2, 3, 4))
getBED(D=dataMZ[[1]], fn=20, ab=3)
}
